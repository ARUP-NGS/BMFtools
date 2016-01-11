#!/usr/bin/env python

import pysam
from enum import Enum
from math import log10
from scipy.stats import combine_pvalues

class UniqueObs(object):
    def __init__(self, BaseCall, FM, PV, FA, MQ, Overlap, Duplex):
        self.BaseCall = BaseCall
        self.FM = FM
        self.FA = FA
        self.PV = PV
        self.MQ = MQ
        self.Overlap = Overlap
        self.Duplex = Duplex

    @classmethod
    def from_name_bin(cls, name_bin):
        if len(name_bin) == 1:
            r = name_bin[0]
            ra = r.alignment
            tagpos = get_tagpos(r)
            FM = ra.opt("FM")
            RV = ra.opt("RV")
            return cls(ra.query_sequence[r.query_position],
                       FM,
                       ra.opt("PV")[tagpos],
                       ra.opt("FA")[tagpos],
                       ra.mapping_quality,
                       False,
                       RV not in [0, FM])
        else:
            r1, r2 = name_bin
            r1a = r1.alignment
            r2a = r2.alignment
            tagpos1 = get_tagpos(r1)
            tagpos2 = get_tagpos(r2)
            bc1 = r1a.query_sequence[r1.query_position]
            bc2 = r2a.query_sequence[r2.query_position]
            if(bc1 == bc2):
                FM = 2 * r1a.opt("FM")
                RV = r1a.opt("RV") + r2a.opt("RV")
                return cls(bc1,
                           2 * r1a.opt("FM"),
                           combine_phreds(r1a.opt("PV")[tagpos1],
                                          r2.Alignment.opt("PV")[tagpos2]),
                           (r1a.opt("FA")[tagpos1] +
                            r2a.opt("FA")[tagpos2]),
                           max(r1a.mapping_quality, r2a.mapping_quality),
                           True, RV not in [0, FM])
            else:
                pv1 = r1a.opt("PV")[tagpos1]
                pv2 = r2a.opt("PV")[tagpos2]
                if(pv1 > pv2):
                    FM = r1a.opt("FM")
                    RV = r1a.opt("RV")
                    return cls(bc1,
                               FM,
                               pv1 - pv2,
                               r1a.opt("FA")[tagpos1],
                               r1a.mapping_quality,
                               False, RV not in [0, FM])
                else:
                    FM = r2a.opt("FM")
                    RV = r2a.opt("RV")
                    return cls(bc2,
                               FM,
                               pv2 - pv1,
                               r2a.opt("FA")[tagpos1],
                               r2a.mapping_quality,
                               False, RV not in [0, FM])
                
                                      

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
    def __init__(self, minFM=0, minPV=0, minFA=0., minMQ=0., minOverlap=0, minCount=0,
                 minDuplex=0):
        self.minFM = minFM
        self.minPV = minPV
        self.minFA = minFA
        self.minMQ = minMQ
        self.minOverlap = minOverlap
        self.minCount = minCount
        self.minDuplex = minDuplex

    def pass_uniobs(self, obs):
        assert isinstance(obs, UniqueObs)
        if (obs.MQ < self.minMQ or
            obs.PV < self.minPV or
            obs.FA < self.minFA or
            obs.FM < self.minFM):
            return False
        return True


def phred2pval(phred):
    return 10 ** (phred * -0.1)


def combine_phreds(phred1, phred2):
    return -10 * log10(combine_pvalues([phred2pval(phred1),
                                        phred2pval(phred2)],
                                       "fisher")[1])


def get_tagpos(pr):
    return (pr.alignment.query_length - 1 - pr.query_position
            if(pr.is_reverse) else pr.query_position)



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
        # 3844 is supp/second/qcfail/duplicate/unmapped
        if i.is_refskip or i.is_del or i.alignment.flag & 3844:
            continue
        if i.query_name in ret:
            ret[i.query_name].append(i)
        else:
            ret[i.query_name] = [i]
    return ret


def get_bc_bins(passing_obs):
    ret = {}
    for i in passing_obs:
        if i.BaseCall in ret:
            ret[i.BaseCall].append(i)
        else:
            ret[i.BaseCall] = [i]
    return ret


class BaseCallSummary(object):
    def __init__(self, bin_tuple, settings, ref_base):
        obs = bin_tuple[1]
        self.is_var = bin_tuple[0] == ref_base
        self.Count = len(obs)
        self.Duplex = sum(i.Duplex for i in obs)
        self.Overlap = sum(i.Overlap for i in obs)
        self.Pass = (self.Duplex > settings.minDuplex and
                     self.Overlap > settings.minOverlap and
                     self.Count > settings.minCount)

def get_bc_bin_summary(bc_bins, vet_settings, ref_base):
        return [BaseCallSummary(bin_tuple, vet_settings, ref_base)
                for bin_tuple in bc_bins.iteritems()]


def get_plp_summary(plp, refstring, vet_settings):
    assert isinstance(vet_settings, vet_set)
    assert isinstance(refstring, str)
    assert isinstance(plp, pysam.calignedsegment.PileupColumn)
    ref_base = refstring[plp.pos]
    observations = map(UniqueObs.from_name_bin,
                       get_name_bins(plp.pileups).itervalues())
    passing = filter(vet_set.pass_unibos, observations)
    bc_bins = get_bc_bins(passing)
    summary = get_bc_bin_summary(bc_bins, vet_settings, ref_base)
    return [i for i in summary if i.is_var and i.Pass]


class Vetter(object):
    def __init__(self, bam, fasta, vcf, settings):
        assert isinstance(bam, pysam.AlignmentFile)
        assert isinstance(fasta, pysam.FastaFile)
        self.settings = settings
        self.bam = bam
        self.fasta = fasta
        self.vcf = vcf

        