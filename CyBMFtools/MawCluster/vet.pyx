#!/usr/bin/env python

import pysam
from enum import Enum
from math import log10
from scipy.stats import combine_pvalues

def get_plp(af=None, contig=None, pos=None):
    """
    af is a pysam.AlignmentFile object.
    contig is a string.
    pos (integer) is 0-based, like in a vcf.
    """
    if pos is None:
        raise ValueError("Need a position argument.")
    if contig is None:
        raise ValueError("Need a contig argument.")
    if af is None:
        raise ValueError("Need an alignmentfile argument.")
    assert isinstance(pos, int)
    assert isinstance(contig, str)
    assert isinstance(af, pysam.AlignmentFile)
    plp_iterator = af.pileup(contig, pos, pos + 1)
    plp_col = plp_iterator.next()
    while plp_col.pos < pos:
        plp_col = plp_iterator.next()
    return plp_col

class filters(Enum):
    minOverlap = 1
    minCount = 2


class vet_set(object):
    def __init__(self, minFM=0, minPV=0, minFA=0., minMQ=0., minOverlap=0, minCount=0):
        self.minFM = minFM
        self.minPV = minPV
        self.minFA = minFA
        self.minMQ = minMQ
        self.minOverlap = minOverlap
        self.minCount = minCount

    def pass_var(self, plpr):
        assert isinstance(plpr, pysam.calignedsegment.PileupRead)
        i = plpr.alignment
        if(i.is_refskip or i.is_del or i.opt("FM") < self.minFM or
           i.mapping_quality < self.minMQ):
            return False
        if self.minPV or self.minFA:
            tagpos = (i.query_length - 1 - i.query_position
                      if(i.is_reverse) else i.query_position)
            if(i.opt("PV")[tagpos] < self.minPV or i.opt("FA")[tagpos] < self.minFA):
                return False
        return True


def phred2pval(phred):
    return 10 ** (phred * -0.1)


def combine_phreds(phred1, phred2):
    return -10 * log10(combine_pvalues([phred2pval(phred1),
                                        phred2pval(phred2)],
                                       "fisher")[1])


class UniqueObs(object):
    def __init__(self, FM, PV, FA, MQ, Overlap):
        self.FM = FM
        self.FA = FA
        self.PV = PV
        self.MQ = MQ
        self.Overlap = Overlap

    @classmethod
    def from_name_bin(cls, name_bin):
        if len(name_bin) == 1:
            r = name_bin[0]
            tagpos = r.alignment.query_length - 1 - r.query_position
            return cls(r.alignment.opt("FM"),
                       r.alignment.opt("PV")[tagpos],
                       r.alignment.opt("FA")[tagpos],
                       r.alignment.mapping_quality,
                       False)
        else:
            r1, r2 = name_bin
            tagpos1 = r1.alignment.query_length - 1 - r1.query_position
            tagpos2 = r2.alignment.query_length - 1 - r2.query_position
            return cls(2 * r1.alignment.opt("FM"),
                       combine_phreds(r1.alignment.opt("PV")[tagpos1],
                                      r2.Alignment.opt("PV")[tagpos2]),
                       (r1.alignment.opt("FA")[tagpos1] +
                        r2.alignment.opt("FA")[tagpos2]),
                       max(r1.alignment.mapping_quality, r2.alignment.mapping_quality),
                       True)
                                      
def freq(iter_obj):
    ret = {}
    for i in iter_obj:
        if i in ret:
            ret[i] += 1
        else:
            ret[i] = 1

def get_name_bins(iter_obj):
    ret = {}
    for i in iter_obj:
        if i.query_name in ret:
            ret[i.query_name].append(i)
        else:
            ret[i.query_name] = [i]
    return ret


def get_plp_summary(plp, refstring, vet_settings):
    flag = 0  # Flag to return which filters need to be applied
    assert isinstance(vet_settings, vet_set)
    assert isinstance(refstring, str)
    assert isinstance(plp, pysam.calignedsegment.PileupColumn)
    ref_base = refstring[plp.pos]
    name_counts = get_name_bins(plp.pileups)
    passing = filter(vet_set.pass_var, plp.pileups)
