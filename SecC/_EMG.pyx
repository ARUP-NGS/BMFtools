# cython: boundscheck=False, wraparound=False
from itertools import groupby
from matplotlib.backends.backend_pdf import PdfPages
from utilBMF.ErrorHandling import ThisIsMadness, ImproperArgumentError
from utilBMF.HTSUtils import pFastqProxy, pFastqFile, getBS, RevCmp, TrimExt
import argparse
import matplotlib.pyplot as plt
import numpy as np
import operator
import pysam
import cython
from array import array
# cimport pysam.calignmentfile
# ctypedef pysam.calignedsegment.AlignedSegment AlignedSegment_t


cdef errorTracker(AlignedSegment_t read,
                  ndarray[uint64_t, ndim=3, mode="c"] readErr,
                  ndarray[uint64_t, ndim=3, mode="c"] readObs,
                  py_array tmpArr):
    cdef int32_t index = 0
    cdef int32_t start, end
    cdef char base
    cdef int8_t phred_index, char_index
    tmpArr = array('B', read.query_sequence)
    index = read.qstart
    for index, base in enumerate(tmpArr[read.qstart:read.qend],
                                 start=read.qstart):
        if base == 78:  # base is '=' or 'N'
            continue
        fprintf(c_stderr,
                "Hey this is base %c with numeric"
                " value %hi and int value %i.",
                base, base, base)
        char_index = NUC_TO_ARRAY(base)
        phred_index = read.query_qualities[index] - 2
        # Note if/elses with same predicate --> switch
        if base != 61:
            readObs[index][phred_index][char_index] += 1
            readErr[index][phred_index][char_index] += 1


def MakeErrorArray(args):

    cdef size_t rLen, index, read_index, qual_index, context_index
    cdef uint64_t fam_range, qcfail, ReadCount, FamSizeFamCount, FM
    cdef double_t r1mean, r2mean
    cdef ndarray[uint64_t, ndim=3, mode="c"] r1error, r1obs
    cdef ndarray[uint64_t, ndim=3, mode="c"] r2error, r2obs
    cdef ndarray[double_t, ndim=3, mode="c"] r1frac, r2frac
    cdef ndarray[double_t, ndim=1, mode="c"] r1cyclemeans, r2cyclemeans
    cdef ndarray[double_t, ndim=1, mode="c"] r1qualmeans, r2qualmeans
    cdef AlignedSegment_t read
    cdef AlignmentFile_t mdBam
    cdef py_array tmpArr

    from cPickle import dump
    from sys import stderr

    if(args.pickle_path is None):
        pickle_path = TrimExt(args.mdBam) + ".errorprofile.pyd"
    else:
        pickle_path = args.pickle_path
    if(args.table_prefix is None):
        table_prefix = TrimExt(args.mdBam)
    else:
        table_prefix = args.table_prefix
    table_prefix += ".out."

    rLen = pysam.AlignmentFile(args.mdBam).next().inferred_length
    tmpArr = array('B')
    c_array.resize(tmpArr, rLen)
    r1error = np.zeros([rLen, 39, 4], dtype=np.uint64)
    r1obs = np.zeros([rLen, 39, 4], dtype=np.uint64)
    r2error = np.zeros([rLen, 39, 4], dtype=np.uint64)
    r2obs = np.zeros([rLen, 39, 4], dtype=np.uint64)

    qcfail = 0
    FamSizeFamCount = 0
    ReadCount = 0

    minFM, maxFM = args.minFM, args.maxFM
    mdBam = pysam.AlignmentFile(args.mdBam, "rb")
    for read in mdBam:
        '''
        First if is equivalent to the following:
            if(read.is_secondary or read.is_supplementary or
               read.is_unmapped or read.is_qcfail):
        Second if is equivalent to not (read.is_proper_pair)
        '''
        if read.flag & 2820 or not (read.flag & 2):
            qcfail += 1
            continue
        FM = read.opt("FM")
        if(FM < minFM or FM > maxFM):
                FamSizeFamCount += 1
                continue
        if read.flag & 64:
            errorTracker(read, r1error, r1obs, tmpArr)
        elif read.flag & 128:
            errorTracker(read, r2error, r2obs, tmpArr)
        else:
            pass
        ReadCount += 1
    stderr.write("Family Size Range: %i-%i\n" % (minFM, maxFM))
    stderr.write("Reads Analyzed: %i\n" % (ReadCount))
    stderr.write("Reads QC Filtered: %i\n" % (qcfail))
    stderr.write("Reads Family Size Filtered: %i\n" % (FamSizeFamCount))
    r1frac = np.divide(r1error, r1obs)  # Now r1frac holds error / obs
    r2frac = np.divide(r2error, r2obs)  # Now r2frac holds error / obs
    r1mean = np.mean(r1frac)
    r2mean = np.mean(r2frac)
    r1cyclemeans = np.mean(np.mean(r1frac, axis=2, dtype=np.float64),
                           axis=1, dtype=np.float64)
    r2cyclemeans = np.mean(np.mean(r2frac, axis=2, dtype=np.float64),
                           axis=1, dtype=np.float64)
    r1qualmeans = np.mean(np.mean(r1frac, axis=0, dtype=np.float64),
                          axis=1, dtype=np.float64)
    r2qualmeans = np.mean(np.mean(r2frac, axis=0, dtype=np.float64),
                          axis=1, dtype=np.float64)
    return {"r1cm": r1cyclemeans, "r2cm": r2cyclemeans,
            "r1qm": r1qualmeans, "r2qm": r2qualmeans,
            "r1f": r1frac, "r2f": r2frac, "r1m": r1mean,
            "r2m": r2mean, "r1e": r1error, "r2e": r2error}


def calculateError(args):

    cdef dict data
    from sys import stderr

    if(args.table_prefix is None):
        table_prefix = TrimExt(args.mdBam)
    else:
        table_prefix = args.table_prefix
    table_prefix += ".out."

    FullTableHandle = open(table_prefix + "full.tsv", "w")
    CycleTableHandle = open(table_prefix + "cycle.tsv", "w")
    PhredTableHandle = open(table_prefix + "phred.tsv", "w")
    ContextTableHandle = open(table_prefix + "context.tsv", "w")

    data = MakeErrorArray(args)
    '''
    CycleTableHandle.write("#Cycle\tRead1 mean error\tRead2 mean error\tread count\n")
    for index in xrange(1, rLen):
        CycleTableHandle.write(
            "%i\t%f\t%f\t%i\n" % (index+1, [index],
                                  r2cyclemeans[index], ReadCount))
    CycleTableHandle.close()
    PhredTableHandle.write("#Quality Score\tr1 mean error\tr2 mean error\n")
    for index in xrange(39):
        PhredTableHandle.write(
            "%i\t%f\t%f\n" % (index+2, r1qualmeans[index],
                              r2qualmeans[index]))
    PhredTableHandle.close()
    ContextTableHandle.write("#Context ID\tRead 1 mean error\tRead 2 mean error\n")
    for index in xrange(16):
        ContextTableHandle.write(
            "%i\t%f\t%f\n" % (index, r1contextmeans[index],
                              r2contextmeans[index]))
    ContextTableHandle.close()
    FullTableHandle.write(
        "#Cycle\tQuality Score\tContext ID\tRead "
        "1 mean error\tRead 2 mean error\n")
    for read_index in xrange(rLen):
        for qual_index in xrange(39):
            for context_index in xrange(16):
                FullTableHandle.write(
                    "%i\t%i\t%i\t%f\t%f\n" % (read_index + 1, qual_index + 2, context_index,
                                              r1frac[read_index][qual_index][context_index],
                                              r2frac[read_index][qual_index][context_index]))
    FullTableHandle.close()
    '''
    return 0


@cython.initializedcheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline bint TEST_ERROR(char character) nogil:
    if character == 61:  # '='
        return 0
    elif character == 78:  # 'N'
        return 0
    else:
        return 1


cpdef errorFinder(AlignedSegment_t read,
                  ndarray[uint64_t, ndim=1] readErr,
                  ndarray[uint64_t, ndim=1] readObs):
    '''
    cdef size_t read_index
    cdef char base
    cdef char * seq
    cdef size_t offset_index
    cdef bint err

    seq = read.query_sequence
    for read_index in xrange(read.qstart, read.qend):
        readObs[read_index] += 1
        err = TEST_ERROR(<char>seq[read_index])
        if not err:
            # case "=" or "N"
            return
        else:
            readErr[read_index] += 1
            return
    '''


@cython.returns(dict)
def cCycleError(args):
    cdef uint64_t rLen, minFM, maxFM, QCFailCount, ReadCount, FamSizeFamCount, index, FM
    cdef double_t r1mean, r2mean
    cdef AlignedSegment_t read
    cdef AlignmentFile_t mdBam
    from sys import stdout, stderr, maxint
    from operator import itemgetter
    cdef ndarray[uint64_t, ndim=1, mode="c"] r1error
    cdef ndarray[uint64_t, ndim=1, mode="c"] r1obs
    cdef ndarray[uint64_t, ndim=1, mode="c"] r2error
    cdef ndarray[uint64_t, ndim=1, mode="c"] r2obs
    rLen = pysam.AlignmentFile(args.mdBam).next().inferred_length
    r1error = np.zeros(rLen, dtype=np.uint64)
    r1obs = np.zeros(rLen, dtype=np.uint64)
    r2error = np.zeros(rLen, dtype=np.uint64)
    r2obs = np.zeros(rLen, dtype=np.uint64)

    if(not args.family_size):  # args.family_size is None
        minFM = 0
        maxFM = maxint
    else:
        minFM, maxFM = itemgetter(0, -1)(args.family_size.split(","))
    QCFailCount = 0
    FamSizeFamCount = 0
    ReadCount = 0
    mdBam = pysam.AlignmentFile(args.mdBam)
    for read in mdBam:
        if read.flag & 2820 or (~read.flag) & 2:
            QCFailCount += 1
            continue
        FM = read.opt("FM")
        if(FM < minFM or FM > maxFM):
                FamSizeFamCount += 1
                continue
        if read.flag & 64:
            errorFinder(read, r1error, r1obs)
        elif read.flag & 128:
            errorFinder(read, r2error, r2obs)
        else:
            pass
        ReadCount += 1
    stdout.write("Family Size Range: %i-%i\n" % (minFM,
                                                 maxFM))
    stdout.write("Reads Analyzed: %i\n" % (ReadCount))
    stdout.write("Reads QC Filtered: %i\n" % (QCFailCount))
    if args.family_size is not None:
        stdout.write("Reads Family Size Filtered: %i\n" % (FamSizeFamCount))
    r1prop = np.divide(r1error.astype(np.double), r1obs)
    r2prop = np.divide(r2error.astype(np.double), r2obs)
    r1mean = np.mean(r1prop)
    r2mean = np.mean(r2prop)
    stdout.write("cycle\tr1\tr2\tread count\n")
    for index in xrange(rLen):
        stdout.write("%i\t%f\t%f\t%i\n" % (index + 1,
                                           r1prop[index],
                                           r2prop[index], ReadCount))
    stdout.write("%i\t%f\t%f\t%i\n" % (maxFM, r1mean, r2mean,
                                       ReadCount))
    if args.cycleheat is not None:
        cycleHeater(r1prop, r2prop, rLen)
    return {"r1e": r1error, "r1o": r1obs, "r2e": r2error,
            "r2o": r2obs, "r1p": r1prop, "r2p": r2prop,
            "r1m": r1mean, "r2m": r2mean}


def cycleHeater(r1prop, r2prop, rLen):
    fig, ax = plt.subplots()
    data = np.array([r1prop, r2prop])
    heatmap = ax.pcolor(data, cmap=plt.cm.Reds)
    colorb = plt.colorbar(heatmap)
    plt.show()
