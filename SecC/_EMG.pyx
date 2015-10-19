# cython: boundscheck=False, wraparound=False
from array import array
from itertools import groupby
from matplotlib.backends.backend_pdf import PdfPages
from utilBMF.ErrorHandling import ThisIsMadness, ImproperArgumentError
from utilBMF.HTSUtils import pFastqProxy, pFastqFile, getBS, RevCmp, TrimExt
from sys import stderr
from libc.string cimport *
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


@cython.boundscheck(False)
@cython.wraparound(False)
cdef errorTracker(AlignedSegment_t mdRead, AlignedSegment_t bamRead,
                  ndarray[uint64_t, ndim=3, mode="c"] readErr,
                  ndarray[uint64_t, ndim=3, mode="c"] readObs,
                  ndarray[double_t, ndim=2, mode="c"] illSum,
                  ndarray[uint64_t, ndim=2, mode="c"] illCount):
    cdef int32_t index = 0
    cdef int32_t start, end
    cdef char base
    cdef int8_t phred_index, char_index
    cdef double_t pVal
    tmpMArr = array('B', mdRead.query_sequence)
    tmpBArr = array('B', bamRead.query_sequence)
    tmpQArr = array('B', bamRead.query_qualities)
    for index, base in enumerate(tmpMArr[mdRead.qstart:mdRead.qend],
                                 start=mdRead.qstart):
        if base == 78:  # base is 'N'
            continue
        char_index = NUC_TO_ARRAY(tmpBArr[index])
        phred_index = tmpQArr[index] - 2
        readObs[index][phred_index][char_index] += 1
        pVal = 10.0**-((phred_index + 2.0)/10.0)
        illSum[index][char_index] += pVal
        illCount[index][char_index] += 1
        if base != 61:
            readErr[index][phred_index][char_index] += 1


cdef pvErrorTracker(AlignedSegment_t read,
                    ndarray[uint64_t, ndim=1, mode="c"] readErr,
                    ndarray[uint64_t, ndim=1, mode="c"] readObs):
    cdef int32_t index = 0
    cdef int32_t start, end
    cdef int32_t phred_index
    cdef double_t pVal
    tmpBArr = array('B', read.query_sequence)
    tmpQArr = read.opt('PV')
    for index, base in enumerate(tmpBArr[read.qstart:read.qend],
                                 start=read.qstart):
        if base == 78:
            continue
        phred_index = tmpQArr[index]
        readObs[phred_index] += 1
        if base != 61:
            readErr[phred_index] += 1


@cython.returns(dict)
def makePVErrorArray(args):
    cdef size_t rLen
    cdef uint64_t qcfail, ReadCount
    cdef ndarray[uint64_t, ndim=1, mode="c"] r1error, r1obs, r2error, r2obs
    cdef ndarray[uint64_t, ndim=1, mode="c"] r1err, r2err
    cdef AlignedSegment_t read
    cdef AlignmentFile_t bam
    r1error = np.zeros(3500, dtype=np.uint64)
    r1obs = np.zeros(3500, dtype=np.uint64)
    r2error = np.zeros(3500, dtype=np.uint64)
    r2obs = np.zeros(3500, dtype=np.uint64)
    bam = pysam.AlignmentFile(args.mdBam)

    qcfail = 0
    ReadCount = 0

    for read in bam:
        if read.flag & 2820 or not (read.flag & 2):
            qcfail+=1
            continue
        if read.flag & 64:
            pvErrorTracker(read, r1error, r2obs)
        elif read.flag & 128:
            pvErrorTracker(read, r2error, r2obs)
        else:
            continue
        ReadCount += 1
    stderr.write("Reads Analyzed: %i\n" % (ReadCount))
    stderr.write("Reads QC Filtered: %i\n" % (qcfail))
    r1err = np.divide(r1error, r1obs)
    r2err = np.divide(r2error, r2obs)
    return {"read1": r1err, "read2": r2err}


@cython.returns(dict)
def MakeErrorArray(args):
    cdef size_t rLen
    cdef uint64_t qcfail, ReadCount
    cdef double_t r1mean, r2mean
    cdef ndarray[uint64_t, ndim=3, mode="c"] r1error, r1obs, r2error, r2obs
    cdef ndarray[double_t, ndim=2, mode="c"] ill1Sum, ill2Sum
    cdef ndarray[uint64_t, ndim=2, mode="c"] ill1Count, ill2Count
    cdef AlignedSegment_t mdRead, bamRead
    cdef AlignmentFile_t mdBam, bam

    bamGen = pysam.AlignmentFile(args.bam).next
    rLen = pysam.AlignmentFile(args.mdBam).next().inferred_length
    r1error = np.zeros([rLen, 39, 4], dtype=np.uint64)
    r1obs = np.zeros([rLen, 39, 4], dtype=np.uint64)
    r2error = np.zeros([rLen, 39, 4], dtype=np.uint64)
    r2obs = np.zeros([rLen, 39, 4], dtype=np.uint64)
    ill1Sum = np.zeros([rLen, 4], dtype=np.float64)
    ill2Sum = np.zeros([rLen, 4], dtype=np.float64)
    ill1Count = np.zeros([rLen, 4], dtype=np.uint64)
    ill2Count = np.zeros([rLen, 4], dtype=np.uint64)

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
            errorTracker(mdRead, bamRead, r1error, r1obs, ill1Sum, ill1Count)
        elif mdRead.flag & 128:
            errorTracker(mdRead, bamRead, r2error, r2obs, ill2Sum, ill2Count)
        else:
            continue
        ReadCount += 1
    stderr.write("Reads Analyzed: %i\n" % (ReadCount))
    stderr.write("Reads QC Filtered: %i\n" % (qcfail))
    return {"read1": {"err": r1error, "obs": r1obs, "illSum": ill1Sum,
            "illCount": ill1Count},
           "read2": {"err": r2error, "obs": r2obs, "illSum": ill2Sum,
            "illCount": ill2Count}}

def format2DOutput(dict data, cystr output):
    """Formats and outputs 2D error array of Cycle and Nucleotide"""
    cdef ndarray[double_t, ndim=2, mode="c"] r1err, r2err
    r1err = np.divide(np.mean(data['read1']['err'].astype(np.float64), axis=1),
                      np.mean(data['read1']['obs'].astype(np.float64), axis=1))
    r2err = np.divide(np.mean(data['read2']['err'].astype(np.float64), axis=1),
                      np.mean(data['read2']['obs'].astype(np.float64), axis=1))
    ill1mean = np.divide(data['read1']['illSum'].astype(np.float64),
                         data['read1']['illCount'].astype(np.float64))
    ill2mean = np.divide(data['read2']['illSum'].astype(np.float64),
                         data['read2']['illCount'].astype(np.float64))
    r1err = -10*np.log10(r1err) - -10*np.log10(ill1mean)
    r2err = -10*np.log10(r2err) - -10*np.log10(ill2mean)
    cdef int cycle, base
    with open(output, 'w') as o:
        for cycle in range(len(r1err)):
            r1bases = []
            r2bases = []
            for base in range(len(r1err[cycle])):
                r1bases.append(str(r1err[cycle][base]))
                r2bases.append(str(r2err[cycle][base]))
            r1bases = ",".join(r1bases)
            r2bases = ",".join(r2bases)
            out = "|".join([r1bases, r2bases])+"\n"
            o.write(out)
        return 0


def genarateCompleteArray(dict data, obsCutoff=1000):
    """creates the complete 3D error array, including filling in
    cycle, qual, base values with inadeuate data with a value derived
    from the differenc between the illumina quality and the
    observed error rates that that nuc/cycle combination"""
    cdef ndarray[double_t, ndim=2, mode="c"] r12d, r122d, ill1mean, ill2mean
    cdef ndarray[double_t, ndim=2, mode="c"] r1offset, r2offset
    cdef ndarray[double_t, ndim=3, mode="c"] r1DataArray, r2DataArray
    cdef ndarray[double_t, ndim=3, mode="c"] r1err, r2err
    r1err = np.divide(data['read1']['err'].astype(np.float64),
                      data['read1']['obs'].astype(np.float64))
    r2err = np.divide(data['read2']['err'].astype(np.float64),
                      data['read2']['obs'].astype(np.float64))
    r12d = np.divide(np.mean(data['read1']['err'].astype(np.float64), axis=1),
                      np.mean(data['read1']['obs'].astype(np.float64), axis=1))
    r22d = np.divide(np.mean(data['read2']['err'].astype(np.float64), axis=1),
                      np.mean(data['read2']['obs'].astype(np.float64), axis=1))
    ill1mean = np.divide(data['read1']['illSum'].astype(np.float64),
                         data['read1']['illCount'].astype(np.float64))
    ill2mean = np.divide(data['read2']['illSum'].astype(np.float64),
                         data['read2']['illCount'].astype(np.float64))
    r1offset = -10*np.log10(r12d) - -10*np.log10(ill2mean)
    r2offset = -10*np.log10(r22d) - -10*np.log10(ill2mean)
    arrShape = [r1err.shape[0], r1err.shape[1], r1err.shape[2]]
    r1DataArray = np.zeros(arrShape, dtype=np.float64)
    r2DataArray = np.zeros(arrShape, dtype=np.float64)
    for cycle in range(len(r1err)):
        for qual in range(len(r1err[cycle])):
            for base in range(len(r1err[cycle][qual])):
                b1 = r1err[cycle][qual][base]
                b2 = r2err[cycle][qual][base]
                b1 = -10*np.log10(b1)
                b2 = -10*np.log10(b2)
                if data['read1']['obs'][cycle][qual][base] < obsCutoff:
                    b1 = qual + 2 + r1offset[cycle][base]
                elif np.isnan(b1) or np.isinf(b1):
                    b1 = qual + 2 + r1offset[cycle][base]
                if data['read2']['obs'][cycle][qual][base] < obsCutoff:
                    b2 = qual + 2 + r2offset[cycle][base]
                if np.isnan(b2) or np.isinf(b2):
                    b2 = qual + 2 + r2offset[cycle][base]
                if b1 <= 0:
                    b1 = 2
                if b2 <= 0:
                    b2 = 2
                r1DataArray[cycle][qual][base] = b1
                r2DataArray[cycle][qual][base] = b2
    return r1DataArray, r2DataArray


def format3DOoutput(ndarray[double_t, ndim=3, mode="c"] r1Array,
                    ndarray[double_t, ndim=3, mode="c"] r2Array,
                    cystr output):
    """Outputs the formatted 3D matrix"""
    with open(output, 'w') as out:
        for ci, cycle in enumerate(r1Array):
            r1quals = []
            r2quals = []
            for qi, qual in enumerate(cycle):
                r1bases = []
                r2bases = []
                for bi, base in enumerate(qual):
                    r1bases.append(str(base))
                    r2bases.append(str(r2Array[ci][qi][bi]))
                r1quals.append(":".join(r1bases))
                r2quals.append(":".join(r2bases))
            r1quals = ",".join(r1quals)
            r2quals = ",".join(r2quals)
            o = "|".join([r1quals, r2quals])+"\n"
            out.write(o)
    return


def calculateErrorArray(args):
    cdef dict data
    from sys import stderr
    if(args.table_prefix is None):
        table_prefix = TrimExt(args.mdBam)
    else:
        table_prefix = args.table_prefix
    table_prefix += ".out"
    data = MakeErrorArray(args)
    if args.dim == "2D":
        format2DOutput(data, table_prefix)
        return 0
    if args.dim == "3D":
        r1dataArray, r2dataArray = genarateCompleteArray(data, args.obsCutoff)
        format3DOoutput(r1dataArray, r2dataArray, table_prefix)
    return 0


def calculatePVErrorArray(args):
    cdef dict data
    if(args.table_prefix is None):
        table_prefix = TrimExt(args.mdBam)
    else:
        table_prefix = args.table_prefix
    table_prefix += ".out"
    data = makePVErrorArray(args)
    plotPVerror(data)
    return 0


def plotPVerror(dict data):
    r1err = data['read1']
    r2err = data['read2']
    fig, ax = plt.subplots()
    plt.plot(np.arange(len(r1err)), r1err, label="read 1")
    plt.plot(np.arange(len(r2err)), r2err, label="read 2")
    plt.legend()
    plt.savefig("wubalubadubdub.png")
    plt.show()

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


def errorFinder():
    pass


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


cdef qualityObservationsArray(cystr inBam, illumina=True):
    cdef py_array qualityCounts = array('L')
    cdef AlignedSegment_t read
    cdef AlignmentFile_t bam
    cdef cystr qual
    cdef uint64_t arrLen
    obs=0
    if illumina == True:
        arrLen = 40
    else:
        arrLen = 4000
    c_array.resize(qualityCounts, arrLen)
    memset(qualityCounts.data.as_voidptr, 0, arrLen*8)
    assert np.sum(qualityCounts) == 0
    bam = pysam.AlignmentFile(inBam, 'rb')
    for read in bam:
        if(illumina == True):
            for i in read.qual:
                qualityCounts[ord(i)-35] += 1
        else:
            for i in read.opt('PV'):
                qualityCounts[i] += 1
    return qualityCounts


def graphQuality(qualArray, name):
    fig, ax = plt.subplots()
    data = np.log10(np.array(qualArray))
    bars = plt.bar(np.arange(len(data)), data)
    plt.savefig(name+".png")
    plt.show()


def qualityRescale(args):
    dmpQualArray = qualityObservationsArray(args.dmpBam, illumina=False)
    inQualArray = qualityObservationsArray(args.inBam, illumina=True)
    if(args.graph == True):
        graphQuality(dmpQualArray, name="dmpQuals")
        graphQuality(inQualArray, name="illuminQuals")
    return 0
