# cython: boundscheck=False, wraparound=False
from __future__ import division
import numpy as np
from utilBMF.HTSUtils import pFastqProxy, pFastqFile, getBS, RevCmp, TrimExt
from itertools import groupby
from utilBMF.ErrorHandling import ThisIsMadness, ImproperArgumentError
from MawCluster.BCFastq import GetDescTagValue
import argparse
import operator
import pysam
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
# cimport pysam.calignmentfile
# ctypedef pysam.calignedsegment.AlignedSegment AlignedSegment_t


cdef errorTracker(AlignedSegment_t read,
                  ndarray[int64_t, ndim=3, mode="c"] readErr,
                  ndarray[int64_t, ndim=3, mode="c"] readObs):
    cdef int32_t index, start, end
    cdef cystr context
    cdef char base
    cdef int8_t phred_index
    cdef int16_t context_index
    for index in xrange(read.qstart, read.qend):
        if(index < 1):
            # Don't sweat it
            continue
        context = read.query_sequence[index - 1: index + 1]
        context_index = CONTEXT_TO_ARRAY_POS(<char *>context)
        if(context_index < 0):
            # Don't sweat it - there was an N in the context
            continue
        phred_index = read.query_qualities[index] - 2
        readObs[index][phred_index][context_index] += 1
        base = read.seq[index]
        # Note if/elses with same predicate --> switch
        if base == 61:
            continue
        elif base == 78:
            continue
        else:
            readErr[index][phred_index][context_index] += 1


def calculateError(args):
    from sys import stderr
    from cPickle import dump
    cdef size_t rLen, index, read_index, qual_index, context_index
    cdef int64_t fam_range, qcfail, rc, fmc, FM
    cdef double_t read1mean, read2mean
    cdef ndarray[int64_t, ndim=3, mode="c"] read1error, read1obs
    cdef ndarray[int64_t, ndim=3, mode="c"] read2error, read2obs
    cdef ndarray[double_t, ndim=3, mode="c"] read1frac, read2frac
    cdef ndarray[double_t, ndim=1, mode="c"] read1cyclemeans, read2cyclemeans
    cdef ndarray[double_t, ndim=1, mode="c"] read1qualmeans, read2qualmeans
    cdef ndarray[double_t, ndim=1, mode="c"] read1contextmeans, read2contextmeans
    cdef AlignedSegment_t read
    cdef AlignmentFile_t mdBam
    if(args.pickle_path is None):
        pickle_path = TrimExt(args.mdBam) + ".errorprofile.pyd"
    else:
        pickle_path = args.pickle_path
    if(args.table_prefix is None):
        table_prefix = TrimExt(args.mdBam)
    else:
        table_prefix = args.table_prefix
    table_prefix += ".out."
    FullTableHandle = open(table_prefix + "full.tsv", "w")
    CycleTableHandle = open(table_prefix + "cycle.tsv", "w")
    PhredTableHandle = open(table_prefix + "phred.tsv", "w")
    ContextTableHandle = open(table_prefix + "context.tsv", "w")
    rLen = pysam.AlignmentFile(args.mdBam).next().inferred_length
    read1error = np.zeros([rLen, 37, 16], dtype=np.int64)
    read1obs = np.zeros([rLen, 37, 16], dtype=np.int64)
    read2error = np.zeros([rLen, 37, 16], dtype=np.int64)
    read2obs = np.zeros([rLen, 37, 16], dtype=np.int64)
    qcfail = 0
    fmc = 0
    rc = 0
    minFM, maxFM = args.minFM, args.maxFM
    mdBam = pysam.AlignmentFile(args.mdBam, "rb")
    for read in mdBam:
        '''
        First if is equivalent to the following:
            if(read.is_secondary or read.is_supplementary or
               read.is_unmapped or read.is_qcfail):
        Second if is equivalent to not (read.is_proper_pair)
        '''
        if read.flag & 2820:
            qcfail += 1
            continue
        elif not (read.flag & 2):
            qcfail += 1
            continue
        FM = read.opt("FM")
        if(FM < minFM or FM > maxFM):
                fmc += 1
                continue
        if read.flag & 64:
            errorTracker(read, read1error, read1obs)
        elif read.flag & 128:
            errorTracker(read, read2error, read2obs)
        else:
            pass
        rc += 1
    stderr.write("Family Size Range: %i-%i\n" % (minFM, maxFM))
    stderr.write("Reads Analyzed: %i\n" % (rc))
    stderr.write("Reads QC Filtered: %i\n" % (qcfail))
    stderr.write("Reads Family Size Filtered: %i\n" % (fmc))
    read1frac = read1error / read1obs
    read2frac = read2error / read2obs
    pickleHandle = open(pickle_path, "wb")
    dump(read1frac, pickleHandle)
    dump(read2frac, pickleHandle)
    pickleHandle.close()
    read1mean = np.mean(read1frac)
    read2mean = np.mean(read2frac)
    read1cyclemeans = np.mean(np.mean(read1frac, axis=2, dtype=np.float64),
                              axis=1, dtype=np.float64)
    read2cyclemeans = np.mean(np.mean(read2frac, axis=2, dtype=np.float64),
                              axis=1, dtype=np.float64)
    read1qualmeans = np.mean(np.mean(read1frac, axis=0, dtype=np.float64),
                              axis=1, dtype=np.float64)
    read2qualmeans = np.mean(np.mean(read2frac, axis=0, dtype=np.float64),
                              axis=1, dtype=np.float64)
    read1contextmeans = np.mean(np.mean(read1frac, axis=0, dtype=np.float64),
                             axis=0, dtype=np.float64)
    read2contextmeans = np.mean(np.mean(read2frac, axis=0, dtype=np.float64),
                             axis=0, dtype=np.float64)
    CycleTableHandle.write("#Cycle\tRead1 mean error\tRead2 mean error\tread count\n")
    for index in xrange(1, rLen):
        CycleTableHandle.write(
            "%i\t%f\t%f\t%i\n" % (index+1, read1cyclemeans[index],
                                  read2cyclemeans[index], rc))
    CycleTableHandle.close()
    PhredTableHandle.write("#Quality Score\tread1 mean error\tread2 mean error\n")
    for index in xrange(37):
        PhredTableHandle.write(
            "%i\t%f\t%f\n" % (index+2, read1qualmeans[index],
                              read2qualmeans[index]))
    PhredTableHandle.close()
    ContextTableHandle.write("#Context ID\tRead 1 mean error\tRead 2 mean error\n")
    for index in xrange(16):
        ContextTableHandle.write(
            "%i\t%f\t%f\n" % (index, read1contextmeans[index],
                              read2contextmeans[index]))
    ContextTableHandle.close()
    FullTableHandle.write(
        "#Cycle\tQuality Score\tContext ID\tRead "
        "1 mean error\tRead 2 mean error\n")
    for read_index in xrange(rLen):
        for qual_index in xrange(37):
            for context_index in xrange(16):
                FullTableHandle.write(
                    "%i\t%i\t%i\t%f\t%f\n" % (read_index + 1, qual_index + 2, context_index,
                                              read1frac[read_index][qual_index][context_index],
                                              read2frac[read_index][qual_index][context_index]))
    FullTableHandle.close()
    return 0

