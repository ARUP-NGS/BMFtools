#!/usr/bin/env python

import pysam
from collections import defaultdict

fail_lines = 0

def var_passes(vrec):
     bp = vrec.info["BMF_VET"]
     return sum(bp) != bp[0]

def generate_passed_lines(vcfpath):
    global fail_lines
    for vrec in pysam.VariantFile(vcfpath):
        if "BMF_VET" in vrec.info:
            if var_passes(vrec):
                yield vrec
            else:
                fail_lines += 1

class VarCounts(object):
    def __init__(self, base, pass_fail, uniobs, qest, uniobs_sum, qest_sum, is_ref=False):
        self.base = base
        self.pass_fail = True if pass_fail else False
        self.uniobs = uniobs
        self.qest = qest
        self.uniobs_freq = uniobs * 1. / uniobs_sum
        self.qest_freq = qest * 1. / qest_sum
        self.is_ref = is_ref

    def __str__(self):
        t = "%c:%s:%i:%i:%f:%f:" % (self.base, self.pass_fail, self.uniobs, self.qest,
                                    self.uniobs_freq, self.qest_freq)
        return t + "ref" if self.is_ref else t + "nonref"

    def is_likely_het(self, minfreq=0.35, maxfreq=0.65):
        return maxfreq > self.qest_freq >= minfreq


class PosVarCounts(object):
    def __init__(self, vrec):
        self.rec = vrec
        self.var = {}
        self.alleles = set()
        self.uniobs_sum = sum(vrec.info['BMF_UNIOBS'])
        self.qest_sum = sum(vrec.info['BMF_QUANT'])
        for allele, val, unic, qest in zip(vrec.alleles,
                                           vrec.info['BMF_VET'],
                                           vrec.info['BMF_UNIOBS'],
                                           vrec.info['BMF_QUANT']):
            if allele == vrec.alleles[0]:
                self.ref = VarCounts(allele, val, unic, qest, self.uniobs_sum, self.qest_sum, True)
            else:
                self.var[allele] = VarCounts(allele, val, unic, qest, self.uniobs_sum, self.qest_sum, False)
            self.alleles.add(allele)

    def get_variant_alleles(self):
        return list(self.var.keys())

    def get_all_alleles(self):
        return list(self.alleles)
                
        


def get_passed_var_freqs(vrec):
    assert var_passes(vrec)
    uniobs_sum = sum(vrec.info['BMF_UNIOBS'])
    qest_sum = sum(vrec.info['BMF_QUANT'])
    ret = {}
    for allele, val, unic, qest in zip(vrec.alleles,
                                       vrec.info['BMF_VET'],
                                       vrec.info['BMF_UNIOBS'],
                                       vrec.info['BMF_QUANT']):
        ret[allele] = VarCounts(allele, val, unic, qest, uniobs_sum, qest_sum,
                                is_ref=allele == vrec.alleles[0])
    print(",".join(map(str, ret.values())))
    return ret


def get_freq_rec_tuples(vrec):
    return vrec, get_passed_var_freqs(vrec)


def generate_passing_var_tuples(vcfpath):
    for vrec in generate_passed_lines(vcfpath):
        yield get_freq_rec_tuples(vrec)
        


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        sys.stderr.write("Usage: %s <vcfpath>.\n" % sys.argv[0])
        sys.exit(1)
    count = 0
    sys.stderr.write("Processing %s...\n" % sys.argv[1])
    for line in generate_passed_lines(sys.argv[1]):
        count += 1
        sys.stdout.write(str(line))
    sys.stderr.write("Wrote %i passing lines, skipped %i failing.\n" % (count, fail_lines))
    sys.exit(0)
