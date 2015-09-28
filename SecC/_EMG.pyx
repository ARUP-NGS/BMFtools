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
import cPickle
# cimport pysam.calignmentfile
# ctypedef pysam.calignedsegment.AlignedSegment AlignedSegment_t


cdef errorTracker(AlignedSegment_t mdRead, AlignedSegment_t bamRead,
                  ndarray[uint64_t, ndim=2, mode="c"] readErr,
                  ndarray[uint64_t, ndim=2, mode="c"] readObs):
    cdef int32_t index = 0
    cdef int32_t start, end
    cdef char base, mbase
    cdef int8_t phred_index, char_index
    tmpMArr = array('B', mdRead.query_sequence)
    tmpBArr = array('B', bamRead.query_sequence)
    for index, base in enumerate(tmpMArr[mdRead.qstart:mdRead.qend],
                                 start=mdRead.qstart):
        if base == 78:  # base is 'N'
            continue
        """
        fprintf(c_stderr,
                "Hey this is base %c with numeric"
                " value %hi and int value %i.",
                base, base, base)
        char_index = NUC_TO_ARRAY(tmpBArr[index])
        phred_index = mdRead.query_qualities[index] - 2
        # Note if/elses with same predicate --> switch
        readObs[index][phred_index][char_index] += 1
        if base != 61:
            readErr[index][phred_index][char_index] += 1
        """
        char_index = NUC_TO_ARRAY(tmpBArr[index])
        readObs[index][char_index] += 1
        if base != 61:
            readErr[index][char_index] += 1


def MakeErrorArray(args):

    cdef size_t rLen, index, read_index, qual_index, context_index
    cdef uint64_t qcfail, ReadCount
    cdef double_t r1mean, r2mean
    cdef ndarray[uint64_t, ndim=2, mode="c"] r1error, r1obs
    cdef ndarray[uint64_t, ndim=2, mode="c"] r2error, r2obs
    cdef ndarray[double_t, ndim=2, mode="c"] r1frac, r2frac
    cdef AlignedSegment_t mdRead, bamRead
    cdef AlignmentFile_t mdBam, bam
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

    bamGen = pysam.AlignmentFile(args.bam).next
    rLen = pysam.AlignmentFile(args.mdBam).next().inferred_length
    tmpArr = array('B')
    c_array.resize(tmpArr, rLen)
    """ preserving 3D arrays here
    r1error = np.zeros([rLen, 39, 4], dtype=np.uint64)
    r1obs = np.zeros([rLen, 39, 4], dtype=np.uint64)
    r2error = np.zeros([rLen, 39, 4], dtype=np.uint64)
    r2obs = np.zeros([rLen, 39, 4], dtype=np.uint64)
    """
    # Moving to 2D arrays here
    r1error = np.zeros([rLen, 4], dtype=np.uint64)
    r1obs = np.zeros([rLen, 4], dtype=np.uint64)
    r2error = np.zeros([rLen, 4], dtype=np.uint64)
    r2obs = np.zeros([rLen, 4], dtype=np.uint64)

    qcfail = 0
    ReadCount = 0

    mdBam = pysam.AlignmentFile(args.mdBam, "rb")
    for mdRead in mdBam:
        bamRead = bamGen()
        '''
        First if is equivalent to the following:
            if(read.is_secondary or read.is_supplementary or
               read.is_unmapped or read.is_qcfail):
        Second if is equivalent to not (read.is_proper_pair)
        '''
        if mdRead.flag & 2820 or not (mdRead.flag & 2):
            qcfail += 1
            continue
        if mdRead.flag & 64:
            errorTracker(mdRead, bamRead, r1error, r1obs)
        elif mdRead.flag & 128:
            errorTracker(mdRead, bamRead, r2error, r2obs)
        else:
            pass
        ReadCount += 1
    stderr.write("Reads Analyzed: %i\n" % (ReadCount))
    stderr.write("Reads QC Filtered: %i\n" % (qcfail))
    #r1frac = np.divide(r1error.astype(np.float64), r1obs.astype(np.float64))
    #r2frac = np.divide(r2error.astype(np.float64), r2obs.astype(np.float64))
    return {'r1o': r1obs, 'r2o': r2obs, 'r1e': r1error, 'r2e': r2error}


def errorArrayFormat2D(dict data, output):
    """formats and outputs error array data in 2D"""
    cdef int cycle, base
    cdef ndarray[uint64_t, ndim=2, mode="c"] r1err, r2err
    r1err = np.divide(data['r1e'].astype(np.float64),
                      data['r1o'].astype(np.float64))
    r2err = np.divide(data['r1e'].astype(np.float64),
                      data['r1o'].astype(np.float64))
    with open(output, 'w') as o:
        for cycle in range(len(r1err)):
            for base in range(len(r2err)):
                b1 = r1err[cycle][base]
                r1bases.append(str(r1err[cycle][base]))
                r2bases.append(str(r2err[cycle][base]))
            r1bases = ",".join(r1bases)
            r2bases = ",".join(r2bases)
            out = "|".join([r1bases, r2bases]) + "\n"
            o.write(out)
    return 0

def errorArrayFormat3D(dict data, output):
    """Formats and outputs error array data"""
    cdef int cycle, qual, base
    with open(output, 'w') as o:
        for cycle in range(len(data["r1f"])):
            r1quals = []
            r2quals = []
            for qual in range(len(data["r1f"][cycle])):
                r1bases = []
                r2bases = []
                for base in range(len(data["r1f"][cycle][qual])):
                    b1 = data['r1f'][cycle][qual][base]
                    b2 = data['r2f'][cycle][qual][base]
                    if np.isnan(b1):
                        b1 = 0.0
                    if np.isnan(b2):
                        b2 = 0.0
                    r1bases.append(str(b1))
                    r2bases.append(str(b2))
                r1quals.append(":".join(r1bases))
                r2quals.append(":".join(r2bases))
            r1quals = ",".join(r1quals)
            r2quals = ",".join(r2quals)
            out = "|".join([r1quals, r2quals])+"\n"
            o.write(out)
    return 0


def calculateErrorArray(args):
    cdef dict data
    from sys import stderr
    if(args.table_prefix is None):
        table_prefix = TrimExt(args.mdBam)
    else:
        table_prefix = args.table_prefix
    table_prefix += ".out"
    """
    FullTableHandle = open(table_prefix + "full.tsv", "w")
    CycleTableHandle = open(table_prefix + "cycle.tsv", "w")
    PhredTableHandle = open(table_prefix + "phred.tsv", "w")
    ContextTableHandle = open(table_prefix + "context.tsv", "w")
    """
    data = MakeErrorArray(args)
    cPickle.dump(data, open(table_prefix+"data.pickle",'w'))
    errorArrayFormat2D(data, table_prefix)
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
