"""
Module for comparing consensus sequences between experiments.

All puns must involve THE LAW.

The Laaaaaaaaaaaaaaaaaaaawwww. LAAAAAAAAAWWWWWWWWWWWWW
"""

import numpy as np
import pysam
import sys

from utilBMF.HTSUtils import TrimExt
from entropy import shannon_entropy as shen


cdef inline cystr GetKeyForRead(AlignedSegment_t read):
    return read.query_name + "_R1" if(read & 64) else read.query_name + "_R2"


cdef inline dict BuildDict(AlignmentFile_t handle, int32_t minFM):
    cdef dict ret
    cdef AlignedSegment_t read
    ret = {GetKeyForRead(read) : (read.opt("FM"), read.query_sequence) for
           read in handle if read.opt("FM") >= minFM}
    return ret


cpdef tuple BuildDictionaries(AlignmentFile_t inHandle1,
                              AlignmentFile_t inHandle2,
                              int32_t minFM):
    """
    Build the hashmaps from barcode to consensus sequences, returning
    both in a tuple.
    """
    cdef tuple ret
    ret = BuildDict(inHandle1, minFM), BuildDict(inHandle2, minFM)
    return ret


cdef inline cystr MakeDiscLine(cystr key, cystr cons1, cystr cons2,
                               int32_t FM1, int32_t FM2, double_t entropy):
    return (key + "\t" + cons1 + "\t" + cons2 + "\t" +
            "%i\t%i\t%f\n" % (FM1, FM2, entropy))


cpdef void WriteDiscToDisk(cystr bam1, cystr bam2, cystr output=None,
                           int32_t minFM=15):
    """
    Finds read families with the same molecular barcode and compares the
    consensus sequences
    :param [cystr/arg] bam1 - Path to first bam.
    :param [cystr/arg] bam2 - Path to second bam.
    :param [cystr/kwarg/None] output - Path to output
    :param [int32_t/kwarg/15] minFM - Minimum family size for
    comparing consensus sequences.
    """
    cdef AlignmentFile_t inHandle1, inHandle2
    cdef dict seqDict1, seqDict2
    cdef cystr key, cons1, cons2
    cdef int32_t FM1, FM2, NumDisc
    inHandle1 = pysam.AlignmentFile(bam1)
    inHandle2 = pysam.AlignmentFile(bam2)
    if output == "stdout":
        outHandle = sys.stdout
    else:
        if output is None:
            output = TrimExt(bam1) + ".discordantpairs.txt"
        outHandle = open(output, "w")
    outHandle.write("#Barcode\tConsensus Sequence #1\tConsensus Sequence #2"
                    "\tFamily Size #1\tFamily Size #2\tShannon Entropy of M"
                    "olecular Barcode\n")
    seqDict1, seqDict2 = BuildDictionaries(inHandle1, inHandle2, minFM)
    for key in sorted(seqDict1.iterkeys()):
        FM1, cons1 = seqDict1[key]
        try:
            FM2, cons2 = seqDict2[key]
        except KeyError:
            sys.stderr.write("Family present in handle 1, but "
                             "not handle 2. Don't sweat it.\n")
            continue
        if cons1 != cons2:
            outHandle.write(MakeDiscLine(key, cons1, cons2,
                                         FM1, FM2, shen(key)))
    outHandle.close()
    return
