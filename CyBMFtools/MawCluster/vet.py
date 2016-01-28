#!/usr/bin/env python

import pysam
from enum import Enum
from math import log10
from scipy.stats import combine_pvalues
from array import array
from sys import float_info
from math import log10 as mlog10

MIN_FLOAT = float_info.min
MAX_PHRED = -10 * mlog10(MIN_FLOAT)

"""
Tools for creating a vetter using metadata from molecular barcodes.

TODO:
  1. Test accessing VCF with pysam.VariantFile
  2. Add the BMF line to the header
  3. Come up with BMF annotation(s).
    1. Pass/Fail
    2. Pass/Fail settings annotations (minCount, minDuplex, &c.)
    3. Estimated quantitation
      1. How about a global expectation for error rate and then ... I don't know.
  4. 
"""

class InfoField(object):
    def __init__(self, id, number, type, description,
                 version=None, source=None):
        self.number = str(number)
        self.type = str(type)
        assert self.type in ["Integer", "Float", "String", "Flag", "Character"]
        self.description = str(description)
        self.id =str(id)

    def __str__(self):
        source_str = ",Source=\"%s\"" % self.source if self.source else ""
        source_str += ",Version=\"%s\"" % self.version if self.version else ""
        return ("##INFO=<ID=%s,Number=%s,"
                "Type=%s,Description=\"%s\"%s>" % (self.id, self.number,
                                                   self.type, self.description,
                                                   source_str))
BMFInfoField = InfoField("BMF", 'A', "String",
                         ("Annotations for pass/fail and quantitation"
                          " estimation using additional molecular "
                          "barcode metadata."))
BMFPass = InfoField("BMF_PASS", 4, "Integer",
                     "Pass/Fail for A/C/G/T variant at position.")

bmf_header_lines = [BMFInfoField, BMFPass]

def get_bmf_info_str():
    return str(BMFInfoField)

class UniqueObs(object):
    def __init__(self, base_call, FM, PV, FA, MQ, overlap, duplex):
        self.base_call = base_call
        self.FM = FM
        self.FA = FA
        self.PV = PV
        self.MQ = MQ
        self.overlap = overlap
        self.duplex = duplex

    @classmethod
    def from_name_bin(cls, name_bin):
        if len(name_bin) == 1:
            r = name_bin[0]
            ra = r.alignment
            rao = ra.opt
            tagpos = get_tagpos(r)
            FM = rao("FM")
            RV = rao("RV")
            return cls(ra.query_sequence[r.query_position],
                       FM,
                       rao("PV")[tagpos],
                       rao("FA")[tagpos],
                       ra.mapping_quality,
                       False,
                       RV not in [0, FM])
        else:
            r1, r2 = name_bin
            r1a = r1.alignment
            r1o = r1a.opt
            r2a = r2.alignment
            r2o = r2a.opt
            tagpos1 = get_tagpos(r1)
            tagpos2 = get_tagpos(r2)
            bc1 = r1a.query_sequence[r1.query_position]
            bc2 = r2a.query_sequence[r2.query_position]
            if(bc1 == bc2):
                FM = 2 * r1o("FM")
                RV = r1o("RV") + r2o("RV")
                return cls(bc1,
                           2 * r1o("FM"),
                           combine_phreds(r1o("PV")[tagpos1],
                                          r2o("PV")[tagpos2]),
                           (r1o("FA")[tagpos1] +
                            r2o("FA")[tagpos2]),
                           max(r1a.mapping_quality, r2a.mapping_quality),
                           True, RV not in [0, FM])
            else:
                pv1 = r1o("PV")[tagpos1]
                pv2 = r2o("PV")[tagpos2]
                if(pv1 > pv2):
                    FM = r1o("FM")
                    RV = r1o("RV")
                    return cls(bc1,
                               FM,
                               pv1 - pv2,
                               r1o("FA")[tagpos1],
                               r1a.mapping_quality,
                               False, RV not in [0, FM])
                else:
                    FM = r2o("FM")
                    RV = r2o("RV")
                    return cls(bc2,
                               FM,
                               pv2 - pv1,
                               r2o("FA")[tagpos1],
                               r2a.mapping_quality,
                               False, RV not in [0, FM])
                

class filters(Enum):
    minOverlap = 1
    minCount = 2


class VetSettings(object):
    def __init__(self, minFM=0, minPV=0, minFA=0., minMQ=0., minOverlap=0, minCount=0,
                 minDuplex=0, fasta=None):
        self.minFM = minFM
        self.minPV = minPV
        self.minFA = minFA
        self.minMQ = minMQ
        self.minOverlap = minOverlap
        self.minCount = minCount
        self.minDuplex = minDuplex
        self.ref = fasta

    def pass_uniobs(self, obs):
        assert isinstance(obs, UniqueObs)
        return (obs.MQ >= self.minMQ and
                obs.PV >= self.minPV and
                obs.FA >= self.minFA and
                obs.FM >= self.minFM)


def phred2pval(phred):
    ret = 10. ** (phred * -0.1)
    return ret


def combine_phreds(phred1, phred2):
    new_pvalue = combine_pvalues([phred2pval(phred1),
                                  phred2pval(phred2)],
                                  "fisher")[1]
    if new_pvalue < MIN_FLOAT:
        return MAX_PHRED
    else:
        return -10 * mlog10(new_pvalue)



def get_tagpos(pr):
    return (pr.alignment.query_length - 1 - pr.query_position
            if(pr.alignment.is_reverse) else pr.query_position)



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
        qname = i.alignment.query_name
        if qname in ret:
            ret[qname].append(i)
        else:
            ret[qname] = [i]
    return ret


def get_bc_bins(passing_obs):
    ret = {}
    for i in passing_obs:
        if i.base_call in ret:
            ret[i.base_call].append(i)
        else:
            ret[i.base_call] = [i]
    return ret


class BaseCallSummary(object):
    """
    BaseCallSummary is a summary of information for a base call
    at a given genomic position.
    """
    def __init__(self, bin_tuple, settings, ref_base):
        obs = bin_tuple[1]
        self.base_call = bin_tuple[0]
        self.count = len(obs)
        dup = 0
        overlap = 0
        for i in obs:
            dup += i.duplex
            overlap += i.overlap
        self.duplex = dup
        self.overlap = overlap
        self.Pass = (self.duplex >= settings.minDuplex and
                     self.overlap >= settings.minOverlap and
                     self.count >= settings.minCount)


def get_bc_bin_summary(bc_bins, vet_settings, ref_base):
        return [BaseCallSummary(bin_tuple, vet_settings, ref_base)
                for bin_tuple in bc_bins.iteritems()]


def get_plp_summary(plp, vet_settings):
    assert isinstance(vet_settings, VetSettings)
    assert isinstance(plp, pysam.calignedsegment.PileupColumn)
    ref_base = vet_settings.ref.fetch(plp.reference_name, plp.pos, plp.pos + 1)
    observations = map(UniqueObs.from_name_bin,
                       get_name_bins(plp.pileups).itervalues())
    print("N observations: %i." % len(observations))
    passing = [obs for obs in observations if vet_settings.pass_uniobs(obs)]
    print("N passing: %i." % len(passing))
    bc_bins = get_bc_bins(passing)
    print("N bins: %i." % len(bc_bins))
    summary = get_bc_bin_summary(bc_bins, vet_settings, ref_base)
    return [i for i in summary if i.Pass]


class Vetter(object):
    def __init__(self, bam, fasta, vcf, settings):
        assert isinstance(bam, pysam.calignmentfile.AlignmentFile)
        assert isinstance(fasta, pysam.cfaidx.FastaFile)
        self.settings = settings
        self.bam = bam
        self.fasta = fasta
        self.vcf = vcf


def add_info_line(header, info_field):
    header.info.add(info_field.id, info_field.number,
                    info_field.type, info_field.description)

def vet_vcf(vf_path, outvf_path, bampath, refpath, outmode="w", **kwargs):
    # cdef pysam.cbcf.VariantFile invf = pysam.cbcf.VariantFile(vf_path)
    # Open variant file handles
    settings = VetSettings(fasta=pysam.FastaFile(refpath), **kwargs)
    invf = pysam.cbcf.VariantFile(vf_path)
    outvf = pysam.cbcf.VariantFile(outvf_path, outmode, header=invf.header)
    bam = pysam.AlignmentFile(bampath, "rb")
    from pysam.cbcf import add_info_array
    # Add header lines
    [add_info_line(outvf.header, line) for line in bmf_header_lines]
    # cdef pysam.cbcf.VariantRecord rec
    ovw = outvf.write
    pass_arr = array('i', [0, 0, 0, 0])
    num_passed = 0
    num_failed = 0
    passing_nucs = set()
    for rec in invf:
        # Goes to pileup at that position
        plp_iterator = bam.pileup(rec.contig, rec.pos - 1)
        plp_col = plp_iterator.next()
        while plp_col.pos < rec.pos:
            plp_col = plp_iterator.next()
        passing_vars = get_plp_summary(plp_col, settings)
        passing_nucs = set([i.base_call for i in passing_vars])
        pass_arr[0] = 'A' in passing_nucs
        pass_arr[1] = 'C' in passing_nucs
        pass_arr[2] = 'G' in passing_nucs
        pass_arr[3] = 'T' in passing_nucs
        print ("Passing variants: '%s'" % ", ".join([i for i in passing_nucs]))
        add_info_array(rec, "BMF_PASS", pass_arr)
        if(sum(pass_arr) != 0):
            "passed"
            num_passed += 1
            ovw(rec)
        else:
            print "Failed"
            num_failed += 1
    invf.close()
    outvf.close()
    return outvf_path
