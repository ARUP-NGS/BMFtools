# cython: c_string_type=str, c_string_encoding=ascii
import abc
from copy import copy as ccopy
from cytoolz import memoize, frequencies as cyfreq
from functools import partial
from itertools import groupby, tee, chain, combinations, product
from MawCluster.Probability import GetCeiling
from numpy import (any as npany, concatenate as nconcatenate, less as nless,
                   max as nmax, mean as nmean, min as nmin)
from operator import iadd as oia, itemgetter as oig, methodcaller as mc
from os import path as ospath
from pysam.calignedsegment import AlignedSegment as pAlignedSegment
from subprocess import check_output, check_call, CalledProcessError
from utilBMF.ErrorHandling import (ThisIsMadness as Tim, FPStr,
                                   FunctionCallException,
                                   IllegalArgumentError, PermissionException,
                                   UnsetRequiredParameter)
from collections import deque
from entropy import shannon_entropy as shen
import copy
import cStringIO
import cython
import logging
import MawCluster
import numconv
import numpy as np
import operator
import os
import pysam
import shlex
import subprocess
import sys
import time
import uuid
import warnings
try:
    from re2 import finditer, compile as regex_compile
except ImportError:
    print("Tried to load re2, a faster version of re. Failed. "
          "Importing regular re.")
    import re
    from re import finditer, compile as regex_compile

oig1 = oig(1)
oig0 = oig(0)
cfi = chain.from_iterable
mcgroup = mc("group", 0)


global __version__

__version__ = "0.1.1"


def l1(x):
    return x[1]


@cython.returns(bint)
def lisnone(x):
    return x[1] is None


@cython.returns(tuple)
def lreverse(x):
    return (x[1], x[0])


@cython.returns(int)
def linsertsize(x):
    return x.insert_size


@cython.returns(int)
def LambdaInsertSize(ReadPair_t x):
    return x.insert_size


def printlog(message, level=logging.INFO):
    message = message.replace(
        "\t", "\\t").replace("'", "\'").replace('"', '\\"')
    Logger = logging.getLogger("Primarylogger")
    if(level == logging.DEBUG):
        # Doesn't print to stderr if set to debug mode.
        Logger.debug(message)
        return
    elif(level == logging.WARNING):
        Logger.warning(message)
        warnings.warn(message)
    elif(level == logging.INFO):
        Logger.info(message)
    else:
        Logger.critical(message)
    sys.stderr.write(message.replace(
        "\t", "\\t").replace(
        "'", "\'").replace('"', '\\"') + "\n")
    return

pl = printlog

nucConverter = numconv.NumConv(4, "ACGT")  # Used to do permutations
nci = nucConverter.int2str
ncs = nucConverter.str2int

nucList = ["A", "C", "G", "T"]


cpdef list permuteNucleotides(long maxn, object nci=nci, int kmerLen=-1):
    """
    nci should be set to a numConv object's int2str method call.
    """
    cdef list tmpList
    cdef size_t tmpInt, strLen
    if(kmerLen < 0):
        return [nci(tmpInt) for tmpInt in xrange(maxn)]
    else:
        tmpList = [nci(tmpInt) for tmpInt in xrange(maxn)]
        for tmpInt in range(maxn):
            strLen = len(tmpList)
            if strLen < kmerLen:
                tmpList[tmpInt] = "A" * (kmerLen - strLen) + tmpList[tmpInt]
        return tmpList


@memoize
@cython.returns(cystr)
def MemoRevCmp(cystr seq):
    return RevCmp(seq)


cpdef cystr RevCmp(cystr seq):
    """
    Reverse complement quickly
    In [48]: %timeit qstr.translate(b)[::-1]
    100000 loops, best of 3: 3.27 us per loop
    # This RevCmp function was the old implicit.
    In [49]: %timeit RevCmp(qstr)
    10000 loops, best of 3: 25.1 us per loop
    """
    return seq.translate(DNA_CODON_TABLE)[::-1]


PysamToChrDict = {}
for i in xrange(22):
    PysamToChrDict[i] = str(i + 1)
PysamToChrDict[-1] = "*"
PysamToChrDict[22] = "X"
PysamToChrDict[23] = "Y"
PysamToChrDict[24] = "MT"
PysamToChrDict[25] = "GL000207.1"
PysamToChrDict[26] = "GL000226.1"
PysamToChrDict[27] = "GL000229.1"
PysamToChrDict[28] = "GL000231.1"
PysamToChrDict[29] = "GL000210.1"
PysamToChrDict[30] = "GL000239.1"
PysamToChrDict[31] = "GL000235.1"
PysamToChrDict[32] = "GL000201.1"
PysamToChrDict[33] = "GL000247.1"
PysamToChrDict[34] = "GL000245.1"
PysamToChrDict[35] = "GL000197.1"
PysamToChrDict[36] = "GL000203.1"
PysamToChrDict[37] = "GL000246.1"
PysamToChrDict[38] = "GL000249.1"
PysamToChrDict[39] = "GL000196.1"
PysamToChrDict[40] = "GL000248.1"
PysamToChrDict[41] = "GL000244.1"
PysamToChrDict[42] = "GL000238.1"
PysamToChrDict[43] = "GL000202.1"
PysamToChrDict[44] = "GL000234.1"
PysamToChrDict[45] = "GL000232.1"
PysamToChrDict[46] = "GL000206.1"
PysamToChrDict[47] = "GL000240.1"
PysamToChrDict[48] = "GL000236.1"
PysamToChrDict[49] = "GL000241.1"
PysamToChrDict[50] = "GL000243.1"
PysamToChrDict[51] = "GL000242.1"
PysamToChrDict[52] = "GL000230.1"
PysamToChrDict[53] = "GL000237.1"
PysamToChrDict[54] = "GL000233.1"
PysamToChrDict[55] = "GL000204.1"
PysamToChrDict[56] = "GL000198.1"
PysamToChrDict[57] = "GL000208.1"
PysamToChrDict[58] = "GL000191.1"
PysamToChrDict[59] = "GL000227.1"
PysamToChrDict[60] = "GL000228.1"
PysamToChrDict[61] = "GL000214.1"
PysamToChrDict[62] = "GL000221.1"
PysamToChrDict[63] = "GL000209.1"
PysamToChrDict[64] = "GL000218.1"
PysamToChrDict[65] = "GL000220.1"
PysamToChrDict[66] = "GL000213.1"
PysamToChrDict[67] = "GL000211.1"
PysamToChrDict[68] = "GL000199.1"
PysamToChrDict[69] = "GL000217.1"
PysamToChrDict[70] = "GL000216.1"
PysamToChrDict[71] = "GL000215.1"
PysamToChrDict[72] = "GL000205.1"
PysamToChrDict[73] = "GL000219.1"
PysamToChrDict[74] = "GL000224.1"
PysamToChrDict[75] = "GL000223.1"
PysamToChrDict[76] = "GL000195.1"
PysamToChrDict[77] = "GL000212.1"
PysamToChrDict[78] = "GL000222.1"
PysamToChrDict[79] = "GL000200.1"
PysamToChrDict[80] = "GL000193.1"
PysamToChrDict[81] = "GL000194.1"
PysamToChrDict[82] = "GL000225.1"
PysamToChrDict[83] = "GL000192.1"
# PysamToChrDict[84] = "gi|9626372|ref|NC_001422.1|"

for i in xrange(22):
    PysamToChrDict[str(i + 1)] = i
PysamToChrDict["*"] = -1
PysamToChrDict["X"] = 22
PysamToChrDict["Y"] = 23
PysamToChrDict["MT"] = 24
PysamToChrDict["GL000207.1"] = 25
PysamToChrDict["GL000226.1"] = 26
PysamToChrDict["GL000229.1"] = 27
PysamToChrDict["GL000231.1"] = 28
PysamToChrDict["GL000210.1"] = 29
PysamToChrDict["GL000239.1"] = 30
PysamToChrDict["GL000235.1"] = 31
PysamToChrDict["GL000201.1"] = 32
PysamToChrDict["GL000247.1"] = 33
PysamToChrDict["GL000245.1"] = 34
PysamToChrDict["GL000197.1"] = 35
PysamToChrDict["GL000203.1"] = 36
PysamToChrDict["GL000246.1"] = 37
PysamToChrDict["GL000249.1"] = 38
PysamToChrDict["GL000196.1"] = 39
PysamToChrDict["GL000248.1"] = 40
PysamToChrDict["GL000244.1"] = 41
PysamToChrDict["GL000238.1"] = 42
PysamToChrDict["GL000202.1"] = 43
PysamToChrDict["GL000234.1"] = 44
PysamToChrDict["GL000232.1"] = 45
PysamToChrDict["GL000206.1"] = 46
PysamToChrDict["GL000240.1"] = 47
PysamToChrDict["GL000236.1"] = 48
PysamToChrDict["GL000241.1"] = 49
PysamToChrDict["GL000243.1"] = 50
PysamToChrDict["GL000242.1"] = 51
PysamToChrDict["GL000230.1"] = 52
PysamToChrDict["GL000237.1"] = 53
PysamToChrDict["GL000233.1"] = 54
PysamToChrDict["GL000204.1"] = 55
PysamToChrDict["GL000198.1"] = 56
PysamToChrDict["GL000208.1"] = 57
PysamToChrDict["GL000191.1"] = 58
PysamToChrDict["GL000227.1"] = 59
PysamToChrDict["GL000228.1"] = 60
PysamToChrDict["GL000214.1"] = 61
PysamToChrDict["GL000221.1"] = 62
PysamToChrDict["GL000209.1"] = 63
PysamToChrDict["GL000218.1"] = 64
PysamToChrDict["GL000220.1"] = 65
PysamToChrDict["GL000213.1"] = 66
PysamToChrDict["GL000211.1"] = 67
PysamToChrDict["GL000199.1"] = 68
PysamToChrDict["GL000217.1"] = 69
PysamToChrDict["GL000216.1"] = 70
PysamToChrDict["GL000215.1"] = 71
PysamToChrDict["GL000205.1"] = 72
PysamToChrDict["GL000219.1"] = 73
PysamToChrDict["GL000224.1"] = 74
PysamToChrDict["GL000223.1"] = 75
PysamToChrDict["GL000195.1"] = 76
PysamToChrDict["GL000212.1"] = 77
PysamToChrDict["GL000222.1"] = 78
PysamToChrDict["GL000200.1"] = 79
PysamToChrDict["GL000193.1"] = 80
PysamToChrDict["GL000194.1"] = 81
PysamToChrDict["GL000225.1"] = 82
PysamToChrDict["GL000192.1"] = 83
# PysamToChrDict["gi|9626372|ref|NC_001422.1|"] = 84


def BuildRefIDDict(cystr bam):
    """
    TODO: Switch over to this vs. the PysamToChrDict usage.
    This way, it defaults to the hard-coded, but automatically
    builds it from the BAM file.
    """
    global ReferenceDict
    ReferenceDict = GetPysamToChrDictFromAlignmentFile(
        pysam.AlignmentFile(bam, "rb"))


@cython.returns(dict)
def getRefIDConversionDict():
    if("ReferenceDict" in globals()):
        return globals()["ReferenceDict"]
    else:
        return PysamToChrDict


@cython.returns(dict)
def GetPysamToChrDictFromAlignmentFile(
        pysam.calignmentfile.AlignmentFile alignmentfileObj):
    """
    Returns a dictionary of pysam reference numbers to contig names.
    """
    return dict(list(enumerate(alignmentfileObj.references)))


@cython.returns(dict)
def GetBidirectionalPysamChrDict(
        pysam.calignmentfile.AlignmentFile alignmentfileObj):
    """
    Returns a dictionary of contig names to pysam reference numbers
    and vice versa - bi-directional.
    """
    refList = list(enumerate(alignmentfileObj.references))
    return dict(map(lreverse, refList) + refList)


cdef class pFastqProxy:
    """
    Python container for pysam.cfaidx.FastqProxy with persistence.
    """

    def __init__(self, cystr comment, cystr quality,
                 cystr sequence, cystr name):
        self.comment = comment
        self.quality = quality
        self.sequence = sequence
        self.name = name

    @classmethod
    def fromFastqProxy(cls, FastqProxy_t FastqProxyObj):
        return cls(FastqProxyObj.comment, FastqProxyObj.quality,
                   FastqProxyObj.sequence, FastqProxyObj.name)

    cdef cystr tostring(self):
        return "@%s %s\n%s\n+\n%s\n" % (self.name, self.comment,
                                        self.sequence, self.quality)

    @cython.returns(cystr)
    def __str__(self):
        return self.tostring()

    cpdef int getFM(self):
        return self.cGetFM()

    cdef int cGetFM(self):
        cdef cystr entry, key, value
        for entry in self.comment.split("|")[1:]:
            key, value = entry.split("=")
            if(key == "FM"):
                return int(value)
        return -1

    cdef cystr cGetBS(self):
        cdef cystr entry, key, value
        for entry in self.comment.split("|")[1:]:
            key, value = entry.split("=")
            if(key == "BS"):
                return value
        return ""

    cpdef cystr getBS(self):
        return self.cGetBS()

    cpdef py_array getQualArray(self):
        return cs_to_ph(self.quality)

    cpdef cystr getSlice(self, int start=0, int end=-1337,
                         cystr addComment=""):
        if(end == -1337):
            return "@%s %s%s\n%s\n+\n%s\n" % (self.name, self.comment,
                                              addComment,
                                              self.sequence[start:],
                                              self.quality[start:])
        return "@%s %s%s\n%s\n+\n%s\n" % (self.name, self.comment,
                                          addComment,
                                          self.sequence[start:end],
                                          self.quality[start:end])


cdef cystr cGetBS(pFastqProxy_t read):
    """
    Portable function for getting the barcode sequence from a marked BMFastq
    """
    cdef cystr entry, key, value
    for entry in read.comment.split("|")[1:]:
        key, value = entry.split("=")
        if(key == "BS"):
            return value
    return ""


cpdef cystr getBS(pFastqProxy_t read):
    """cpdef wrapper of cGetBS
    """
    return cGetBS(read)


@cython.returns(cystr)
def FastqProxyToStr(FastqProxy_t fqPrx):
    """
    Just makes a string from a FastqProxy object.
    """
    return "@%s %s\n%s\n+\n%s\n" % (fqPrx.name, fqPrx.comment,
                                    fqPrx.sequence, fqPrx.quality)


@cython.returns(cystr)
def SliceFastqProxy(FastqProxy_t fqPrx,
                    int firstBase=0,
                    int lastBase=-1337,
                    cystr addString=""):
    if(lastBase == -1337):
        return "@%s %s%s\n%s\n+\n%s\n" % (fqPrx.name, fqPrx.comment,
                                          addString,
                                          fqPrx.sequence[firstBase:],
                                          fqPrx.quality[firstBase:])
    return "@%s %s%s\n%s\n+\n%s\n" % (fqPrx.name, fqPrx.comment,
                                      addString,
                                      fqPrx.sequence[firstBase:lastBase],
                                      fqPrx.quality[firstBase:lastBase])


def FacePalm(string, art=FPStr):
    raise Tim("WHAT YOU SAY")


@cython.returns(bint)
def is_read_softclipped(read):
    """
    Simply returns whether or not a read is soft-clipped
    """
    if(read.cigarstring is None):
        return False
    if("S" in read.cigarstring):
        return True
    return False


@cython.returns(bint)
def ReadPairIsDuplex(readPair, minShare="default"):
    """
    If minShare is an integer, require that many nucleotides
    overlapping to count it as duplex.
    Defaults to sharing at least half.
    """
    cdef int minLen
    if(readPair.read1_contig != readPair.read2_contig):
        return False
    if(isinstance(minShare, int)):
        minLen = minShare
    elif(isinstance(minShare, float)):
        minLen = int(minShare * readPair.read1.query_length)
    elif(minShare == "default"):
        minLen = readPair.read1.query_length // 2
    else:
        raise Tim("minShare parameter required. Integer for an absolute "
                  "number of bases overlapped required, float for a frac"
                  "tion of read length.")
    return sum([x == 2 for x in
                cyfreq(readPair.read1.get_reference_positions() +
                       readPair.read2.get_reference_positions()).values()]
               ) >= minLen


def BwaswCall(fq1, fq2, ref="default", outBAM="default"):
    if(ref == "default"):
        raise Tim("ref required to call bwasw.")
    if(outBAM == "default"):
        outBAM = ".".join(fq1.split(".")[:-1]) + ".bam"
    cStr = "bwa bwasw %s %s %s | samtools view -Sbh - > %s" % (ref, fq1, fq2,
                                                               outBAM)
    pl("About to call bwasw. Command string: %s" % cStr)
    check_call(cStr, shell=True)
    return outBAM


def BedtoolsBam2Fq(BAM, outfq1="default", outfq2="default"):
    """
    Converts a BAM to 2 fastq files.
    """
    if(outfq1 == "default"):
        outfq1 = ".".join(BAM.split(".")[:-1] + ["bam2fq.R1.fastq"])
    if(outfq2 == "default"):
        outfq2 = ".".join(BAM.split(".")[:-1] + ["bam2fq.R2.fastq"])
    commandString = "bedtools bamtofastq -i %s -fq %s -fq2 %s" % (
        BAM, outfq1, outfq2)
    check_call(shlex.split(commandString))
    return outfq1, outfq2


def align_bwa_aln(cystr R1, cystr R2, cystr ref=None,
                  cystr opts="", cystr outBAM=None,
                  bint addRG=True):
    """
    Aligns a set of paired-end reads using bwa aln. Defaults to 4 threads.
    In order to make BAMs compatible with both GATK and pysam,

    """
    if(ref is None):
        raise Tim("Reference file index required for alignment!")
    if(opts == ""):
        opts = "-n 3 -t 4"
    if(outBAM is None):
        outBAM = '.'.join(R1.split('.')[0:-1]) + ".aln.bam"
    outSAM = outBAM.replace("bam", "sam")
    # Note: ID has to be "bwa" so that it passes the SAM validation that
    # the hsdjdk has, which is required for either GATK realignment or
    # ABRA realignment.
    R1Sai = R1 + ".tmp.sai"
    R2Sai = R2 + ".tmp.sai"
    alnStr1 = ("bwa aln " + opts + " " + " ".join([ref, R1]) +
               " > " + R1Sai)
    check_call(alnStr1, shell=True)
    alnStr2 = ("bwa aln " + opts + " " + " ".join([ref, R2]) +
               " > " + R2Sai)
    check_call(alnStr2, shell=True)
    sampeBaseStr = "bwa sampe " + " ".join([ref, R1Sai, R2Sai, R1, R2])
    if(addRG):
        sampeBaseStr += (" | sed 's/^@PG/@RG\tID:default\tPL:ILLUMINA\tPU:def"
                         "ault\tLB:default\tSM:default\tCN:default\n@PG/'")
    # sampeStr = sampeBaseStr + " | samtools view -Sbh - > %s" % outBAM
    sampeStr = sampeBaseStr + " | samtools view -Sbh - > %s" % outBAM
    printlog("bwa aln string: {}".format(sampeStr))
    check_call(sampeStr.replace("\t", "\\t").replace("\n", "\\n"), shell=True)
    check_call(["rm", R1Sai, R2Sai])
    return outBAM


def align_bwa_mem(R1, R2, ref="default", opts="", outBAM="default",
                  path="default",
                  bint addCO=True, bint fqCO=True):
    """
    Aligns a set of paired-end
    reads to a reference
    with provided options using bwa mem.
    Defaults to 4 threads, silent alignment, listing
    supplementary alignments, and
    writing each reads' alignment,
    regardless of mapping quality.
    In addition, adds an RG header line for "default",
    primarily for compatibility with GATK/Picard.
    :param cystr R1 - Path to Fq 1
    :param cystr R2 - Path to Fq 2
    :param cystr ref - Path to reference
    :param cystr opts - Options to pass to bwa
    :param bint addCO - Whether or not to use the -C option.
    :param bint fqCO - True if the fastq comment section has
    had CO: prepended to it.
    """
    if(path == "default"):
        path = "bwa"
    if(opts == ""):
        opts = '-t 4 -v 1 -Y -T 0'
    if(outBAM == "default"):
        outBAM = ".".join(R1.split(".")[0:-1]) + ".mem.bam"
    if(ref == "default"):
        raise Tim("Reference file index required for alignment!")
    opt_concat = ' '.join(opts.split())
    baseString = "%s mem %s %s %s %s " % (path, opt_concat, ref, R1, R2)
    if(addCO):
        baseString = baseString.replace("%s mem" % path, "%s mem -C" % path)
        sedString = (" | sed -r -e 's/\t~#!#~[1-4]:[A-Z]:[0-9]+:[AGCNT]+\|/\t"
                     "RG:Z:default\tCO:Z:|/' -e 's/^@PG/@RG\tID:default\tPL:"
                     "ILLUMINA\tPU:default\tLB:default\tSM:default\tCN:defaul"
                     "t\n@PG/'")
        baseString += sedString
    if(path == "default"):
        command_str = baseString + " | samtools view -Sbh - > %s" % outBAM
    else:
        command_str = "%s%s | samtools view -Sbh - > %s" % (path,
                                                            baseString[3:],
                                                            outBAM)
    # command_list = command_str.split(' ')
    printlog("bwa mem command string with RG/CO additions"
             ": %s" % command_str)
    check_call(command_str.replace("\n", "\\n").replace("\t", "\\t"),
               shell=True)
    printlog("bwa mem aligned output is: %s" %outBAM)
    return outBAM


@cython.returns(cystr)
def PipeAlignTag(R1, R2, ref="default",
                 outBAM="default", path="default",
                 bint coorsort=True, bint u=False,
                 cystr sortMem="6G", cystr opts=None,
                 bint dry_run=False, bint sam=False):
    """
    :param R1 - [cystr/arg] - path to input fastq for read 1
    :param R2 - [cystr/arg] - path to input fastq for read 2
    :param ref [cystr/kwarg/"default"] path to reference index base
    :param outBAM - [cystr/kwarg/"default"] - path to output bam.
    Set to 'stdout' to emit to stdout.
    :param path - [cystr/kwarg/"default"] - absolute path to bwa executable.
    :param coorsort [bint/kwarg/True] - whether or not to coordinate sort
    :param u [bint/kwarg/False] - emit uncompressed bam.
    Override default bwa path (bwa) if necessary.
    :param sortMem - [cystr/kwarg/"6G"] - sort memory limit for samtools
    :param opts - [cystr/kwarg/"-t 4 -v 1 -Y -T 0"] - optional arguments
    to provide to bwa for alignment.
    :param dry_run - [bint/kwarg/False] - flag to return the command string
    rather than calling it.
    :returns: [cystr] - path to outBAM if writing to file, "stdout" if
    emitting to stdout.
    """
    if(opts is None):
        opts = "-t 4 -v 1 -Y -T 0"
    if(path == "default"):
        path = "bwa"
    if(outBAM == "default"):
        outBAM = ".".join(R1.split(".")[0:-1]) + ".mem.bam"
    if(ref == "default"):
        raise Tim("Reference file index required for alignment!")
    PBTflag = 6 if(u) else 2
    uuidvar = str(uuid.uuid4().get_hex().upper()[0:8])
    opt_concat = ' '.join(opts.split())
    cStr = "%s mem -C %s %s %s %s " % (path, opt_concat, ref, R1, R2)
    sedString = (" | sed -r -e 's/\t~#!#~[1-4]:[A-Z]:[0-9]+:[AGCNT]*\|/\t"
                 "RG:Z:default\t/' -e 's/^@PG/@RG\tID:default\tPL:"
                 "ILLUMINA\tPU:default\tLB:default\tSM:default\tCN:defaul"
                 "t\n@PG/' -e 's/FP=/FP:i:/' -e 's/\|BS=/\tBS:Z:/' -e "
                 "'s/\|FM=/\tFM:i:/' -e 's/\|ND=/\tND:i:/' -e 's/\|FA=/\t"
                 "FA:B:i,/' -e 's/\|PV=/\tPV:B:i,/'")
    cStr += sedString
    if(sam is False):
        if(coorsort):
            compStr = " -l 0 " if(u) else ""
            cStr += " | samtools sort -m %s -O bam -T %s %s -" % (sortMem,
                                                                  uuidvar,
                                                                  compStr)
            if(outBAM != "stdout" and outBAM != "-"):
                cStr += " -o %s" % outBAM
        else:
            cStr += (" | samtools view -Sbhu - " if(
                u) else " | samtools view -Sbh -")
            if(outBAM != "stdout" and outBAM != "-"):
                cStr += " > %s" % outBAM
    else:
        cStr += " > %s" % outBAM.replace(".bam", ".sam")
    pl("Command string for ambitious pipe call: %s" % cStr.replace(
        "\t", "\\t").replace("\n", "\\n"))
    if(dry_run):
        return cStr.replace("\n", "\\n").replace("\t", "\\t")
    else:
        check_call(cStr.replace("\t", "\\t").replace("\n", "\\n"),
                   shell=True, executable="/bin/bash")
        if(coorsort is True and outBAM != "stdout" and outBAM != "-"):
            check_call(shlex.split("samtools index %s" % outBAM))
        return outBAM


def align_bwa_mem_se(reads, ref, opts, outBAM):
    """Aligns a set of reads to a reference
    with provided options. Defaults to
    4 threads, silent alignment, listing
    supplementary alignments, and
    writing each reads' alignment,
    regardless of mapping quality.
    """
    if(opts == ""):
        opts = '-t 4 -v 1 -Y -T 0'
    opt_concat = ' '.join(opts.split())
    command_str = ('bwa mem {} {} {}'.format(opt_concat, ref, reads) +
                   " | samtools view -Sbh - > {}".format(outBAM))
    # command_list = command_str.split(' ')
    printlog(command_str)
    check_call(command_str, shell=True)
    return outBAM, command_str


def align_snap(R1, R2, ref, opts, outBAM):
    opt_concat = " ".join(opts.split())
    command_str = "snap paired {} {} {} -o {} {}".format(
        ref,
        R1,
        R2,
        outBAM,
        opt_concat)
    printlog(command_str)
    subprocess.check_call(shlex.split(command_str), shell=False)
    return(command_str)


def CustomRefBowtiePaired(mergedFq,
                          ref,
                          output="default",
                          barLen="default",
                          bowtiePath="bowtie",
                          mismatchLimit="default"):
    if(output == "default"):
        output = mergedFq.split('.')[0] + '.mergingFamilies.sam'
    if(barLen == "default"):
        raise Tim("Barcode length must be set. Abort mission!")
    if(bowtiePath == "bowtie"):
        printlog("Defaulting to bowtie for path to executable.")
    elif("2" in bowtiePath.split("/")[-1]):
        raise Tim("Do not use bowtie2!")
    command_list = [
        "bowtie",
        "--threads",
        "4",
        "-S",
        "-n",
        str(mismatchLimit),
        "-l",
        str(barLen),
        "--norc",
        ref,
        mergedFq,
        output]
    subprocess.check_call(command_list)
    return output


# This function is to handle StopIterations with a little elegance
def has_elements(iterable):
    from itertools import tee
    iterable, any_check = tee(iterable)
    try:
        any_check.next()
        return True, iterable
    except StopIteration:
        return False, iterable


def IntervalOverlapsBed(queryInterval, bedIntervals, bedDist=0):
    """
    Requires bedIntervals in the form of the output of ParseBed
    Returns True or False as to whether or not an overlap exists
    for this interval in the bed file.
    Now expanded to create an allowance for distance from bed regions.
    Default behavior has an bedDist of 0, equivalent to the original
    function.
    """
    if(queryInterval[1] > queryInterval[2]):
        newInt = copy.copy(queryInterval)
        newInt[1] = copy.copy(queryInterval[2])
        newInt[2] = copy.copy(queryInterval[1])
        queryInterval = newInt
    for interval in bedIntervals:
        if(queryInterval[0] == interval[0]):
            if(queryInterval[1] > interval[2] - 1 + bedDist or
               queryInterval[2] < interval[1] - 1 + bedDist):
                continue
            else:
                return True
        else:
            continue
    return False


def ReadWithinDistOfBedInterval(samRecord, bedLine="default", dist=70):
    """
    Checks to see if a samRecord is contained in a bedfile.
    bedLine must be a list, where list[0] is a string and
    list[1] and list[2] are integers. ParseBed returns a list of such objects.
    """
    try:
        contig = PysamToChrDict[samRecord.reference_id]
    except KeyError:
        # Read most likely unmapped.
        return False
    if(contig == bedLine[0]):
        if((samRecord.reference_start > bedLine[2] - 1 + dist) or
           (samRecord.reference_end < bedLine[1] - 1 - dist)):
            return False
        else:
            return True
    else:
        return False


def ReadOverlapsBed(samRecord, bedRef="default"):
    """
    Checks to see if a samRecord is contained in a bedfile.
    bedRef must be a tab-delimited list of lists, where
    line[1] and line[2] are integers. ParseBed returns such an object.
    """
    # if(isinstance(bedRef, str)):
    #     bedRef = ParseBed(bedRef)
    assert isinstance(bedRef[0][0], str) and isinstance(bedRef[0][1], int)
    for line in bedRef:
        """
        try:
            assert isinstance(line, list)
        except AssertionError:
            print(repr(line))
            raise Tim("OMGZ")
        """
        try:
            contig = PysamToChrDict[samRecord.reference_id]
        except KeyError:
            # Read most likely unmapped.
            return False
        if(contig == line[0]):
            if((samRecord.reference_start > line[2] - 1) or
               (samRecord.reference_end < line[1] - 1)):
                continue
            else:
                # print("Read {} was contained in bed file".format(
                #       samRecord.query_name))
                return True
        else:
            continue
    # print("Read {} which was not contained in bed file".format(
    #       samRecord.query_name))
    return False


def VCFLineContainedInBed(VCFLineObject, bedRef="default"):
    """
    Checks to see if a VCF Line is contained in a bedfile.
    bedRef must be a tab-delimited list of lists, where
    line[1] and line[2] are integers. ParseBed returns such an object.
    """
    try:
        assert(isinstance(VCFLineObject, MawCluster.BCVCF.VCFRecord))
    except AssertionError:
        raise Tim("VCFLineContainedInBed requires a VCFRecord object!")
    for line in bedRef:
        if(VCFLineObject.CHROM == line[0]):
            assert isinstance(line[1], int)
            assert isinstance(line[2], int)
            assert isinstance(line[0], str)
            assert isinstance(VCFLineObject.POS, int)
            if(VCFLineObject.POS > line[2] or VCFLineObject.POS <= line[1]):
                continue
            else:
                return True
        else:
            continue
    return False


def PosContainedInBed(contig, pos, bedRef="default"):
    """
    Checks to see if a position is contained in a bedfile.
    0-based
    """
    if(isinstance(bedRef, str)):
        pl("Bed file being parsed for each call of PosContainedInBed: "
           "WARNING!")
        bedRef = ParseBed(bedRef)
    for line in bedRef:
        if(contig == line[0]):
            if(pos >= line[2] or pos <= line[1]):
                continue
            else:
                return True
        else:
            continue
    # print("Read {} which was not contained in bed file".format(
    #       samRecord.query_name))
    return False


def indexBowtie(fasta):
    subprocess.check_call('bowtie-build {0} {0}'.format(fasta), shell=True)
    return


def BedtoolsGenomeCov(inBAM, ref="default", outfile="default"):
    if(ref == "default"):
        raise Tim("A reference file path must be provided!")
    if(outfile == "default"):
        outfile = inBAM[0:-3] + ".doc.txt"
    outfileHandle = open(outfile, "w")
    subprocess.check_call(
        (shlex.split("bedtools genomecov -ibam {}".format(inBAM) +
                     " -dz -g {}").format(ref)), stdout=outfileHandle)
    outfileHandle.close()
    return outfile


def BedtoolsBamToBed(inBAM, outbed="default", ref="default"):
    if(ref == "default"):
        raise Tim("A reference file path must be provided!")
    if(outbed == "default"):
        outbed = inBAM[0:-4] + ".doc.bed"
    outfile = BedtoolsGenomeCov(inBAM, ref=ref)
    OutbedAppendingList = []
    lastPos = 0
    outbedHandle = open(outbed, "w")
    for line in [l.strip().split(
                 '\t') for l in open(outfile, "r").readlines()]:
        if(len(OutbedAppendingList) == 0):
            OutbedAppendingList = [line[0], line[1], "unset"]
            lastPos = int(line[1])
        if(int(line[1]) - lastPos == 1):
            lastPos += 1
        else:
            OutbedAppendingList[2] = int(line[1]) + 1
            outbedHandle.write("\t".join(OutbedAppendingList) + "\n")
            OutbedAppendingList = []
    outbedHandle.close()
    return outbed


def CoorSortAndIndexBam(inBAM, prefix="MetasyntacticVar",
                        outBAM="default",
                        uuid="true", memStr="4G",
                        threads="4", delete=False):
    '''
    Uses samtools >= 1.0.0 to coordinate sort and index a bam file.
    If uuid is either a boolean true or is a string containing true,
    then a random string is generated for the output
    '''
    if("true" in str(uuid).lower()):
        import uuid
        prefix += str(uuid.uuid4().get_hex().upper()[0:8])
    if(outBAM == "default"):
        outBAM = '.'.join(inBAM.split('.')[0:-1]) + '.CoorSort.bam'
    CommandStr = ("samtools sort -T {} -O bam -o {}".format(prefix, outBAM) +
                  " -@ {} {} -m {}".format(threads, inBAM, memStr))
    printlog("About to call sort command: {}".format(CommandStr))
    subprocess.check_call(shlex.split(CommandStr))
    printlog("Now indexing.")
    subprocess.check_call(shlex.split("samtools index {}".format(outBAM)))
    if(delete):
        subprocess.check_call(["rm", inBAM])
    printlog("Output BAM is: " + outBAM)
    return outBAM


def NameSort(inBAM, outBAM="default", prefix="MetasyntacticVar",
             uuid="true", threads="4", memStr="4G"):
    # If uuid is either a boolean true or is a string containing true,
    # then a random string is generated for the output
    if(str(uuid).lower() == "true"):
        import uuid
        prefix += str(uuid.uuid4().get_hex().upper()[0:8])
    if(outBAM == "default"):
        outBAM = '.'.join(inBAM.split('.')[0:-1]) + '.NameSort.bam'
    CommandStr = ("samtools sort -T {} -O bam -o {}".format(prefix, outBAM) +
                  " -@ {} -m {} -n {}".format(threads, memStr, inBAM))
    printlog("About to call sort command: {}".format(CommandStr))
    subprocess.check_call(shlex.split(CommandStr))
    printlog("Namesort successful, sorted bam available at: {}".format(outBAM))
    return outBAM


@cython.locals(sortAndIndex=bint)
def NameSortAndFixMate(inBAM, outBAM="default", prefix="MetasyntacticVar",
                       uuid="true", threads="4", memStr="4G",
                       sortAndIndex=False, deleteIntermediate=True):
    """
    If uuid is either a boolean true or is a string containing true,
    then a random string is generated for the output
    """
    tempOutBam = NameSort(inBAM)
    if(outBAM == "default"):
        outBAM = ".".join(inBAM.split(".")[0:-1] + ["ns", "fxm", "bam"])
    subprocess.check_call(shlex.split(
        "samtools fixmate %s %s -O bam" % (tempOutBam, outBAM)))
    if(sortAndIndex):
        return CoorSortAndIndexBam(inBAM, prefix="NSFxM8CS",
                                   delete=deleteIntermediate)
    return outBAM


def mergeBamPicardOld(
        samList, memoryStr="-XmX16",
        MergeJar="/mounts/bin/picard-tools/MergeSamFiles.jar",
        outBAM="default"):
    if(outBAM == "default"):
        outBAM = '.'.join(samList[0].split('.')[0:-1]) + '.merged.bam'
    cStr = ("java -jar " + MergeJar + " " + memoryStr + " I=" +
            " I=".join(samList) + " O=" + outBAM + " MSD=True " +
            "AS=True SO=coordinate"
            )
    printlog("About to merge bams. Command string: " + cStr)
    subprocess.check_call(shlex.split(cStr))
    return outBAM


def samtoolsMergeBam(bamlist, outBAM="default", NameSort=True):
    if(outBAM == "default"):
        outBAM = ".".join(bamlist[0].split(".")[:-1]) + ".merged.bam"
        if(len(outBAM) > 100):
            outBAM = bamlist[0].split(".")[0] + ".merged.bam"
            pl("Shortened output bam name.", level=logging.DEBUG)
    if(NameSort):
        cStr = "samtools merge -n %s %s" % (outBAM, " ".join(bamlist))
    else:
        cStr = "samtools merge %s %s" % (outBAM, " ".join(bamlist))
    pl("Merge bams command: %s" % (cStr))
    check_output(cStr)
    return outBAM


cdef class ReadPair:

    """
    Holds both bam record objects in a pair.
    Currently, one read unmapped and one read soft-clipped are
    both marked as soft-clipped reads.
    """

    def __init__(self, AlignedSegment_t read1,
                 AlignedSegment_t read2):
        self.read1 = read1
        self.read2 = read2
        try:
            self.SVTags = read1.opt("SV").split(',')
        except KeyError:
            self.SVTags = None
        self.insert_size = abs(read1.tlen)
        if(read1.is_unmapped):
            self.read1_is_unmapped = True
            self.read1_soft_clipped = True
        else:
            self.read1_is_unmapped = False
            if("S" in read1.cigarstring):
                self.read1_soft_clipped = True
            else:
                self.read1_soft_clipped = False
        if(read2.is_unmapped):
            self.read2_is_unmapped = True
            self.read2_soft_clipped = True
        else:
            self.read2_is_unmapped = False
            if("S" in read2.cigarstring):
                self.read2_soft_clipped = True
            else:
                self.read2_soft_clipped = False
        if(self.read1_is_unmapped):
            self.read1_contig = "*"
        else:
            self.read1_contig = PysamToChrDict[read1.reference_id]
        if(self.read2_is_unmapped):
            self.read2_contig = "*"
        else:
            self.read2_contig = PysamToChrDict[read2.reference_id]
        self.SameContig = (read1.reference_id == read2.reference_id)
        self.ContigString = ",".join(sorted([self.read1_contig,
                                             self.read2_contig]))
        self.SameStrand = (self.SameContig and
                           (read1.is_reverse == read2.is_reverse))

    @cython.returns(int)
    def NumOverlappingBed(self, list bedLines=[]):
        try:
            assert isinstance(bedLines[0], str) and isinstance(
                bedLines[1], int)
        except AssertionError:
            raise Tim("Sorry, bedLines must be in ParseBed format.")
        if(self.read1_is_unmapped is False):
            self.read1_in_bed = ReadOverlapsBed(self.read1, bedLines)
        else:
            self.read1_in_bed = False
        if(self.read2_is_unmapped is False):
            self.read2_in_bed = ReadOverlapsBed(self.read2, bedLines)
        else:
            self.read2_in_bed = False
        return sum([self.read2_in_bed, self.read1_in_bed])

    @cython.returns(list)
    def getReads(self):
        return [self.read1, self.read2]


@cython.returns(list)
def GetOverlappingBases(ReadPair_t pair):
    """
    Returns the bases of the reference to which both reads in the pair
    are aligned.
    """
    if(pair.SameContig):
        return [i[0] for i in
                cyfreq([pair.read1.aligned_pairs +
                        pair.read2.aligned_pairs]).iteritems() if
                i[1] == 1 and i[0] != None]
    return []


@cython.returns(dict)
def AlignPairDict(AlignedSegment_t read):
    return {x: y for y, x in read.aligned_pairs}


@cython.returns(AlignedSegment_t)
def CollapseReadPair(ReadPair_t pair, bint BMFTags=True,
                     int minQualDiff=3):
    """
    minQualDiff is the minimum difference between the quality
    scores in the case of disagreement.
    """
    cdef AlignedSegment_t read1, read2, newread
    cdef int i
    if(not pair.SameContig or pair.SameStrand):
        return None  # Nothing to collapse!
    read1, read2 = pair.getReads()
    overlap = GetOverlappingBases(pair)
    r1matchdict = AlignPairDict(pair.read1)
    r2matchdict = AlignPairDict(pair.read2)
    r1positions = [r1matchdict[i] for i in overlap]
    r2positions = [r2matchdict[i] for i in overlap]
    r1baseTuples = [(i, read1.seq[i], read1.query_qualities[i]) for
                    i in r1positions]
    r2baseTuples = [(i, read2.seq[i], read2.query_qualities[i]) for
                    i in r2positions]
    CollapsedSeq = ""
    CollapsedNewQuals = []
    newread = pysam.AlignedSegment()
    newread.qname = read1.qname + "Combined"
    newread.is_reverse = read1.is_reverse
    newread.reference_id = read1.reference_id
    newread.pos = min([read1.pos, read2.pos])
    newread.mapq = max([read1.mapq, read2.mapq])
    newread.tlen = read1.tlen
    newread.set_tags(read1.get_tags())
    read1First = (read1.pos > read2.pos)
    if(not BMFTags):
        for r1tup, r2tup in zip(r1baseTuples, r2baseTuples):
            if(r1tup[1] == r2tup[1]):
                CollapsedSeq += r1tup[1]
            else:
                if(r2tup[2] - r1tup[2] > minQualDiff):
                    CollapsedSeq += r2tup[1]
                    CollapsedNewQuals.append(r2tup[2] - r1tup[2])
                elif(r1tup[2] - r2tup[2] > minQualDiff):
                    CollapsedSeq += r1tup[1]
                    CollapsedNewQuals.append(r1tup[2] - r2tup[2])
                else:
                    CollapsedSeq += "N"
                    CollapsedNewQuals.append(0)
    if(read1First):
        raise Tim("Sorry, I just have other things to finish now.")


def CollapseR1R2(AlignedSegment_t R1,
                 AlignedSegment_t R2):
    return CollapseReadPair(ReadPair(R1, R2))


cdef class pPileupRead:
    """
    Python container for the PileupRead proxy in pysam
    """

    cpdef opt(self, cystr arg):
        return self.alignment.opt(arg)

    def __init__(self, pysam.calignedsegment.PileupRead PileupRead):
        cdef py_array BQs
        self.alignment = PileupRead.alignment
        self.indel = PileupRead.indel
        self.level = PileupRead.level
        self.query_position = PileupRead.query_position
        self.name = self.alignment.qname
        self.BaseCall = self.alignment.seq[self.query_position]
        BQs = self.alignment.opt("PV")
        self.BQ = BQs[self.query_position]
        self.MBQ = nmax(BQs)
        BQs = self.alignment.opt("FA")
        self.FA = BQs[self.query_position]
        self.FM = self.alignment.opt("FM")
        self.MQ = self.alignment.mapping_quality


cdef class PileupReadPair:

    """
    Holds both bam record objects in a pair of pileup reads.
    Currently, one read unmapped and one read soft-clipped are
    both marked as soft-clipped reads.
    Accepts a list of length two as input.
    """

    def MarkReads(self):
        for read in self.RP.getReads():
            read.set_tag("DP", self.RP.discordanceString, "Z")

    def __cinit__(self, tuple readlist):
        cdef pPileupRead_t read1
        cdef pPileupRead_t read2
        read1, read2 = readlist[0], readlist[1]
        try:
            assert len(readlist) == 2
        except AssertionError:
            pl("repr(readlist): %s" % repr(readlist))
            raise Tim(
                "readlist must be of length two to make a PileupReadPair!")
        self.RP = ReadPair(read1.alignment, read2.alignment)
        self.read1 = read1
        self.read2 = read2
        self.discordant = (read1.BaseCall != read2.BaseCall)
        self.name = read1.alignment.query_name
        if(self.discordant):
            if(read1.alignment.is_reverse):
                self.discordanceString = (self.RP.read1_contig + "," +
                                          str(self.read1.alignment.pos -
                                              self.read1.query_position))
            else:
                self.discordanceString = (self.RP.read1_contig + "," +
                                          str(self.read1.alignment.pos +
                                              self.read1.query_position))
        else:
            self.discordanceString = ""
        self.MarkReads()


def GetReadPair(inHandle):
    """
    Simply contains both pairs of reads in an object
    """
    read1 = inHandle.next()
    read2 = inHandle.next()
    try:
        assert read1.query_name == read2.query_name
    except AssertionError:
        raise Tim("These two reads have "
                  "different query names. Abort!")
    return ReadPair(read1, read2)


cdef bint cReadsOverlap(AlignedSegment_t read1,
                        AlignedSegment_t read2):
    # Same strand or different contigs.
    if(read1.is_unmapped or read2.is_unmapped or
       read1.reference_id != read2.reference_id or
       read1.is_reverse == read2.is_reverse or
       read1.reference_start > read2.reference_end or
       read1.reference_end < read2.reference_start):
        return False
    if(read1.aend < read2.pos or read2.aend < read1.pos):
        return False
    return True


def ReadPairsInsertSizeWithinDistance(ReadPair1, ReadPair2, distance=300):
    if(abs(ReadPair1.insert_size - ReadPair2.insert_size <= distance)):
        return True
    return False


def ReadPairsOverlap(ReadPair1, ReadPair2):
    if(ReadsOverlap(ReadPair1.read1, ReadPair2.read1) or
       ReadsOverlap(ReadPair1.read2, ReadPair2.read2) or
       ReadsOverlap(ReadPair1.read1, ReadPair2.read2) or
       ReadsOverlap(ReadPair1.read2, ReadPair2.read1)):
        return True
    else:
        return False


def ReadPairsWithinDistance(ReadPair1, ReadPair2, distance=300):
    if(sum(
        [abs(ReadPair1.read1.pos - ReadPair2.read1.pos) < distance,
         abs(ReadPair1.read1.pos - ReadPair2.read2.pos) < distance,
         abs(ReadPair1.read2.pos - ReadPair2.read1.pos) < distance,
         abs(ReadPair1.read2.pos - ReadPair2.read2.pos) < distance]) > 0):
        return True
    else:
        return False


def ReadPairPassesMinQ(Pair, minMQ=0, minBQ=0):
    assert isinstance(Pair, ReadPair)
    if(Pair.read1.mapq < minMQ or Pair.read2.mapq < minMQ):
        return False
    if(npany(nless(Pair.read1.query_qualities, minBQ)) or
       npany(nless(Pair.read2.query_qualities, minBQ))):
        return False
    else:
        return True


def LoadReadsFromFile(inBAM, SVTag="default", minMQ=0,
                      minFamSize="default"):
    RecordsArray = []
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    while True:
        try:
            RecordsArray.append(inHandle.next())
        except StopIteration:
            break
    RecordsArray = [rec for rec in RecordsArray if rec.mapq >= minMQ]
    if(SVTag != "default"):
        for tag in SVTag.split(','):
            RecordsArray = [rec for rec in RecordsArray if tag
                            in rec.opt("SV")]
    if(minFamSize != "default"):
        try:
            minFamSize = int(minFamSize)
        except ValueError:
            raise Tim("Minimum family size must be castable to int!")
        RecordsArray = [rec for rec in RecordsArray if
                        rec.opt("FM") >= minFamSize]
    inHandle.close()
    return RecordsArray


def LoadReadPairsFromFile(inBAM, SVTag="default",
                          minMQ=0, minBQ=0,
                          LambdaInsertSize=LambdaInsertSize):
    """
    Loads all pairs of reads from a name-sorted paired-end
    bam file into ReadPair objects. If SVTag is specified,
    then check that all entries in SVTag.split(",") are in
    the tags
    """
    RecordsArray = []
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    tags = SVTag.split(',')
    print("Tags: {}".format(repr(tags)))
    if(SVTag != "default"):
        while True:
            try:
                read1 = inHandle.next()
                read2 = inHandle.next()
                WorkingReadPair = ReadPair(read1, read2)
                if(WorkingReadPair.read1.mapq >= minMQ and
                   WorkingReadPair.read2.mapq >= minMQ and
                   sum([tag in WorkingReadPair.SVTags
                        for tag in tags]) == len(tags)):
                    RecordsArray.append(WorkingReadPair)
                else:
                    pass
            except StopIteration:
                # print("Stopping iterations...")
                break
    else:
        while True:
            try:
                RecordsArray.append(GetReadPair(inHandle))
            except StopIteration:
                break
    if("LI" in tags):
        return sorted(RecordsArray, key=LambdaInsertSize)
    else:
        return RecordsArray


cpdef bint WritePairToHandle(
        ReadPair_t pair,
        pysam.calignmentfile.AlignmentFile handle=None):
    """
    Writes a pair to a file handle.
    """
    try:
        handle.write(ReadPair.read1)
        handle.write(ReadPair.read2)
        return True
    except Exception:
        return False


@cython.returns(cystr)
def BedListToStr(list bedlist):
    return bedlist[0] + ":%s-%s" % (bedlist[1], bedlist[2])


@cython.returns(list)
def ParseBed(cystr bedfile):
    """
    Parses a bedfile in, leaving a list of length 3.
    bed[0] is a string (contig), and bed[1:] are all
    integers.
    """
    bed = [line.strip().split(
    )[0:3] for line in open(
        bedfile, "r").readlines() if line[0] != "#"]
    for line in bed:
        line[1] = int(line[1])
        line[2] = int(line[2])
    return bed


@cython.returns(bint)
@cython.locals(input_str=cystr)
def to_bool(input_str):
    return (input_str.lower() == "true")

TypeConversionDict = {"s": str, "i": int, "f": float, "b": to_bool}


@cython.locals(lst=list, typechar=cystr,
               TypeConversionDict=dict)
def parseTuple(lst, TypeConversionDict=TypeConversionDict):
    assert(len(lst) == 2)
    try:
        typechar = lst[1][0]
    except IndexError:
        return lst[0]  # Is a string
    return TypeConversionDict[typechar](lst[0])


@cython.locals(path=cystr, parsedLines=list)
@cython.returns(dict)
def parseSketchConfig(path):
    """
    Parses in a file into a dictionary of key value pairs.

    Note: config style is key|value|typechar, where typechar
    is 'b' for bool, 'f' for float, 's' for string, and 'i' for int.
    Anything after a # character is ignored.
    """
    parsedLines = [l.strip().split("#")[0].split("|") for l in
                   open(path, "r").readlines()
                   if l[0] != "#"]
    # Note that the key is mangled to make the key match up with
    # argparse's name
    return {line[0].replace(" ", "_"): parseTuple([line[1], line[2]]) for
            line in parsedLines}


@cython.returns(dict)
def parseConfig(cystr string):
    """
    Parses in a file into a dictionary of key value pairs.
    Key is line.strip().split("=")[0].
    Value is line.strip().split("=")[1].
    Any further values are ignored.
    New with BMFTools v0.0.5.2 (or so?): # comment the rest of a line out.
    """
    parsedLines = [l.strip().split("#")[0] for l in
                   open(string, "r").readlines()
                   if l[0] != "#"]
    return {line.split("=")[0].strip(): line.split("=")[1].strip() for
            line in parsedLines}


@cython.returns(dict)
def ReadListToCovCounter(reads, int minClustDepth=3,
                         int minPileupLen=10):
    """
    Makes a Counter object of positions covered by a set of reads.
    Only safe at this point for intrachromosomal rearrangements!
    """
    return cyfreq(reduce(lambda x, y: x + y,
                         [r.get_reference_positions() for r in reads]))


@cython.returns(dict)
def ReadPairListToCovCounter(list ReadPairList, int minClustDepth=5,
                             int minPileupLen=10):
    """
    Makes a Counter object of positions covered by a set of read pairs.
    Only safe at this point for intrachromosomal rearrangements!
    We discount the "duplex" positions because we want to look for pileups of
    read pairs (ultimately, for supporting a structural variant).
    """
    cdef dict PosCounts
    posList = []
    posListDuplex = []
    for pair in ReadPairList:
        R1Pos = pair.read1.get_reference_positions()
        R2Pos = pair.read2.get_reference_positions()
        oia(posList, R1Pos)
        oia(posList, R2Pos)
        oia(posListDuplex, [pos for pos in R1Pos if pos in R2Pos])
    PosCounts = cyfreq(posList)
    PosDuplexCounts = cyfreq(posListDuplex)
    # decrement the counts for each position to account for
    # both reads in a pair mapping to the same location.
    for key in PosDuplexCounts.iterkeys():
        PosCounts[key] -= PosDuplexCounts[key]
    PosCounts = dict([i for i in PosCounts.iteritems()
                      if i[1] >= minClustDepth])
    return PosCounts


class SoftClippedSeq:
    def __init__(self, seq, contig=None, is_reverse=None):
        assert isinstance(seq, str)
        if isinstance(contig, str) is False:
            raise Tim("Soft-clipped seq is ambiguous without a contig.")
        if isinstance(is_reverse, bool) is False:
            raise Tim("SoftClippedSeq needs to know to "
                      "which strand it is mapped.")
        self.seq = seq
        self.contig = contig
        self.is_reverse = is_reverse

    def RCChangeStrand(self):
        self.seq = RevCmp(self.seq)
        self.is_reverse = not self.is_reverse


def SplitSCRead(AlignedSegment_t read):
    try:
        assert "S" in read.cigarstring
    except AssertionError:
        raise Tim("You can't split a read by soft-clipping "
                  "if it's not soft-clipped!")
    SoftClips = []
    clippedStart = 0
    clippedEnd = read.query_length
    if(read.cigar[0][0] == 4):
        SoftClips.append(SoftClippedSeq(read.str[0:read.cigar[0][1]],
                                        contig=PysamToChrDict[
                                            read.reference_id],
                                        is_reverse=read.is_reverse))
        clippedStart = read.cigar[0][1]
    if(read.cigar[-1][0] == 4):
        SoftClips.append(SoftClippedSeq(read.str[read.cigar[-1][1]:],
                                        contig=PysamToChrDict[
                                            read.reference_id],
                                        is_reverse=read.is_reverse))
        clippedEnd -= read.cigar[-1][1]
    newQuals = read.seq[clippedStart:clippedEnd]
    read.seq = read.seq[clippedStart:clippedEnd]
    read.query_qualities = newQuals
    return read, SoftClips


class Interval:
    def __init__(self, intervalList, meanDOR=0):
        assert isinstance(intervalList[0], str)
        assert isinstance(intervalList[1], int)
        self.interval = intervalList
        self.meanDOR = meanDOR
        self.length = self.interval[2] - self.interval[1]
        self.contig = self.interval[0]
        self.start = self.interval[1]
        self.end = self.interval[2]


@cython.returns(int)
def LambdaSub(x):
    return x[0] - x[1]


@cython.returns(list)
def CreateIntervalsFromCounter(dict CounterObj, int minPileupLen=0,
                               cystr contig="default",
                               bedIntervals="default",
                               int mergeDist=0,
                               int minClustDepth=5, lix=LambdaSub):
    """
    From a dictionary object containing the sum of the output of
    get_reference_positions for a list of AlignedSegment objects, it creates a
    list of continuous intervals, calculates the mean coverage for each, and
    returns two lists, the first containing the 0-based open intervals and
    second containing the mean coverage of that interval.
    bedIntervals must be in ParseBed output format.
    """
    cdef list IntervalList
    cdef list MergedInts
    cdef list MeanCovList
    cdef list posList
    cdef list interval
    IntervalList = []
    MeanCovList = []
    if(contig == "default"):
        raise Tim("contig required for this function!")
    for k, g in groupby(
            enumerate(sorted(CounterObj.iterkeys())), lambda x: x[0] - x[1]):
        posList = map(oig1, g)
        if(posList[0] < posList[-1]):
            interval = [contig, posList[0], posList[-1] + 1]
        else:
            interval = [contig, posList[-1], posList[0] + 1]
        if(interval[2] - interval[1] < minPileupLen):
            print("Interval too short. Length: {}".format(interval[2] -
                                                          interval[1]))
            continue
        IntervalList.append(interval)
        MeanCovList.append(nmean([CounterObj[key] for key in posList if
                                 len(posList) >= minPileupLen]))
    # Now merge, in case some are entirely adjacent
    if(len(IntervalList) == 0):
        return []
    print("Now attempting to merge any adjacent intervals. Number: {}".format(
        len(IntervalList)))
    IntervalList = sorted(IntervalList, key=oig1)
    workingIntval = copy.copy(IntervalList[0])
    MergedInts = []
    for intval in IntervalList:
        if(workingIntval == intval):
            continue
        if(intval[1] - workingIntval[2] < 2 + mergeDist):
            MergedInts.append([workingIntval[0], workingIntval[1], intval[2]])
        else:
            MergedInts.append(workingIntval)
            workingIntval = copy.copy(intval)
    MergedInts.append(intval)
    print("MergedInts={}".format(MergedInts))
    return MergedInts

ph2chrDict = {i: chr(i + 33) if i < 94 else "~" for i in xrange(100000)}
# Pre-computes
chr2ph = {i: ord(i) - 33 for i in [ph2chrDict[i] for i in range(94)]}
chr2phStr = {x: str(y) for x, y in chr2ph.iteritems()}
int2Str = {i: str(i) for i in xrange(10000)}

"""
@cython.returns(np.int64_t)
def chr2ph(x):
    \"""
    Converts a character to its corresponding phred integer representation
    \"""
    return ord(x) - 33
"""


@memoize
def ph2chr(x):
    """
    Converts a phred score to a fastq-encodable character.
    """
    return chr(x + 33) if x <= 93 else "~"


@cython.returns(list)
def GetDeletedCoordinates(AlignedSegment_t read):
    """
    Returns a list of integers of genomic coordinates for deleted bases.
    """
    cdef int k
    apList = read.get_aligned_pairs()
    k = 0
    # Remove any soft-clipped bases at the front.
    while True:
        if(apList[k][0] is None):
            k += 1
        else:
            break
    apList = apList[k:]
    k = len(apList) - 1
    # Remove any soft-clipped bases at the end.
    while True:
        if(apList[k][0] is None):
            k -= 1
        else:
            break
    apList = apList[0:k + 1]
    return sorted([i[1] for i in apList if i[0] is None])


@cython.returns(list)
def GetInsertedNucleotides(AlignedSegment_t read):
    """
    """
    cdef int start, end
    cdef list apList
    cdef tuple i
    apList = read.get_aligned_pairs()
    # Remove any soft-clipped bases at the front.
    genPos = []
    readPos = []
    start = 0
    end = len(apList)
    for k, g in groupby(apList, l1):
        genPos.append(list(g))
        readPos.append(k)
    if(readPos[0] is None):
        start = len(genPos[0])
    if(readPos[-1] is None):
        end = len(apList) - len(genPos[-1])
    apList = apList[start:end]
    return sorted([i[0] for i in apList if i[1] is None])


@cython.returns(list)
def GetInsertedStrs(AlignedSegment_t read):
    """
    Returns a list of tuples for strings, along with preceding reference base
    and successive reference base position. We can then directly compare these
    tuples between read 1 and read 2 if they cover the same positions.
    Used to determine whether or not DSI is appropriate.
    """
    cdef dict readPosToAlignedPosDict
    cdef int l
    cdef list PrecedingBase, SuccessiveBase, stringList, set
    positions = GetInsertedNucleotides(read)
    stringList = []
    PrecedingBase = []
    SuccessiveBase = []
    readPosToAlignedPosDict = dict(read.get_aligned_pairs())

    for k, g in groupby(enumerate(positions), lisnone):
            set = map(oig1, g)
            if(set[0] > set[-1]):
                PrecedingBase.append(readPosToAlignedPosDict[set[-1] - 1])
                SuccessiveBase.append(readPosToAlignedPosDict[set[0] + 1])
            else:
                PrecedingBase.append(readPosToAlignedPosDict[set[0] - 1])
                SuccessiveBase.append(readPosToAlignedPosDict[set[-1] + 1])
            stringList.append("".join([read.seq[l] for l in set]))
    return zip(stringList, PrecedingBase, SuccessiveBase)


@cython.returns(np.longdouble_t)
def FractionSoftClipped(AlignedSegment_t read):
    """
    Returns the fraction of a read aligned.
    """
    if(read.cigarstring is None):
        return 0.
    return FractionSoftClippedCigar(read.cigar)


@cython.returns(np.longdouble_t)
def FractionSoftClippedCigar(list cigar):
    """
    Returns the fraction softclipped directly from tuple
    """
    cdef tuple i
    return 1. * sum(i[1] for i in cigar if i[0] == 4) / sum([i[1] for i
                                                             in cigar])


@cython.returns(np.longdouble_t)
def FractionAlignedCigar(list cigar):
    """
    Returns the fraction aligned directly from tuple
    """
    cdef tuple i
    return 1. * sum(i[1] for i in cigar if i[0] == 0) / sum([i[1] for i
                                                             in cigar])


@cython.returns(np.longdouble_t)
def FractionAligned(AlignedSegment_t read):
    """
    Returns the fraction of a read aligned.
    """
    if(read.cigarstring is None):
        return 0.
    return FractionAlignedCigar(read.cigar)


def AddReadGroupsPicard(inBAM, RG="default", SM="default",
                        PL="ILLUMINA", CN="default", picardpath="default",
                        outBAM="default", ID="default", LB="default",
                        PU="default"):
    if(picardpath == "default"):
        raise Tim("picardpath required to call PicardTools!")
    if(outBAM == "default"):
        outBAM = ".".join(inBAM.split(".")[:-1] + ["addRG", "bam"])
    commandStr = ("java -jar %s AddOrReplaceReadGroups I=" % picardpath +
                  "%s O=%s VALIDATION_STRINGENCY=SILENT " % (inBAM, outBAM) +
                  " CN=%s PL=%s SM=%s ID=%s LB=%s PU=%s" % (CN, PL,
                                                            SM, ID, LB,
                                                            PU))
    printlog("AddReadGroupsPicard commandStr: %s" % commandStr)
    subprocess.check_call(shlex.split(commandStr))
    return outBAM


@cython.locals(outliers_fraction=np.longdouble_t,
               contamination=np.longdouble_t,
               window=int)
def BuildEEModels(f1, f2, outliers_fraction=0.1, contamination=0.005,
                  window=20):
    from sklearn.covariance import EllipticEnvelope
    cdef ndarray[np.longdouble_t, ndim=1] GAFreqNP = f1
    cdef ndarray[np.longdouble_t, ndim=1] CTFreqNP = f2
    cdef ndarray[np.longdouble_t, ndim=1] FreqArray
    FreqArray = nconcatenate(GAFreqNP, CTFreqNP)
    ee1 = EllipticEnvelope(contamination=contamination, assume_centered=False)
    ee2 = EllipticEnvelope(contamination=contamination, assume_centered=False)
    GAClassifier = ee1.fit(GAFreqNP)
    CTClassifier = ee2.fit(CTFreqNP)
    pass


def PlotNormalFit(array, outfile="default"):
    import matplotlib as mpl
    mpl.use('Agg')      # With this line = figure disappears
    import matplotlib.pyplot as plt
    import matplotlib.mlab as mlab
    array = array.reshape(-1, 1)
    mu, sigma = nmean(array), cyStdFlt(array)
    n, bins, patches = plt.hist(array, 50, normed=1, facecolor='green',
                                alpha=0.75)
    y = mlab.normpdf(bins, mu, sigma)
    l = plt.plot(bins, y, 'r--', linewidth=1)
    plt.xlabel('A->C/G->/T frequency')
    plt.ylabel('Proportion of sites')
    plt.title(r'$\mathrm{Histogram\ of\ deamination substitutions:}\ '
              r'\mu=%s,\ \sigma=%s$' % (mu, sigma))
    plt.axis([nmin(array), nmax(array), 0., nmax(n)])
    plt.grid(True)
    plt.savefig(outfile + ".png")
    return outfile


@cython.returns(tuple)
def CalculateFamStats(inFq):
    """
    Calculates mean family size, families of size 1, mean family size for
    families of size > 1, and number of singletons.
    """
    inHandle = pysam.FastqFile(inFq)
    cdef FastqProxy_t read
    cdef int famS
    cdef int sumFam
    cdef int numFam
    cdef int numSing
    cdef int sumAll
    cdef np.longdouble_t meanFamAll
    cdef np.longdouble_t meanRealFam
    numSing = 0
    sumFam = 0
    numFam = 0
    sumSing = 0
    sumAll = 0

    for read in inHandle:
        famS = int(read.comment.split(" ")[3].split("=")[1])
        if(famS > 1):
            numFam += 1
            sumFam += famS
        else:
            numSing += 1
        sumAll += famS
    meanFamAll = sumAll / (1. * numFam + numSing)
    meanRealFam = sumFam / (1. * numFam)
    printlog("Number of reads in familes size >=2 %s" % sumFam)
    printlog("Number of reads in familes, all: %s" % sumAll)
    printlog("Numbers of singleton families: %s." % (numSing))
    printlog("Numbers of families of size >= 2: %s." % (numFam))
    printlog("Mean families size, all: %s." % (meanFamAll))
    printlog("Mean families size, FM >= 2: %s." % (meanRealFam))
    return numSing, numFam, meanFamAll, meanRealFam


@cython.locals(n=int)
@cython.returns(list)
def bitfield(n):
    """
    Parses a bitwise flag into an array of 0s and 1s.
    No need - use the & or | tools for working with bitwise flags.
    """
    return [1 if digit == '1' else 0 for digit in bin(n)[2:]]


@cython.returns(cystr)
def ASToFastqSingle(AlignedSegment_t read):
    """
    Makes a string containing a single fastq record from
    one read.
    """
    if read.is_read1:
        if(read.is_reverse):
            return ("@" + read.query_name + " 1\n" + read.seq +
                    "\n+\n" +
                    RevCmp("".join([ph2chrDict[i] for
                                    i in read.query_qualities])))
        return ("@" + read.query_name + " 1\n" + read.seq +
                "\n+\n" +
                "".join([ph2chrDict[i] for i in read.query_qualities]))
    if read.is_read2:
        if(read.is_reverse):
            return ("@" + read.query_name + " 2\n" + read.seq +
                    "\n+\n" +
                    RevCmp("".join([ph2chrDict[i] for
                                    i in read.query_qualities])))
        return ("@" + read.query_name + " 2\n" + read.seq +
                "\n+\n" +
                "".join([ph2chrDict[i] for i in read.query_qualities]))
    else:
        return ("@" + read.query_name + "\n" + read.seq +
                "\n+\n" +
                "".join([ph2chrDict[i] for i in read.query_qualities]))


@cython.locals(alignmentfileObj=pysam.calignmentfile.AlignmentFile)
@cython.returns(cystr)
def ASToFastqPaired(AlignedSegment_t read,
                    alignmentfileObj):
    """
    Works for coordinate-sorted and indexed BAM files, but throws
    an error if the mate is unmapped.
    """
    FastqStr1 = ASToFastqSingle(read)
    FastqStr2 = ASToFastqSingle(alignmentfileObj.mate(read))
    if(read.is_read1):
        return FastqStr1 + "\n" + FastqStr2
    return FastqStr2 + "\n" + FastqStr1


@cython.returns(pysam.calignmentfile.AlignmentFile)
def SWRealignAS(AlignedSegment_t read,
                pysam.calignmentfile.AlignmentFile alignmentfileObj,
                cystr extraOpts="",
                cystr ref="default", float minAF=0.5):
    """
    Passes the sequence and qualities of the provided read to bwa aln with
    an AlignedSegment as input. Updates the AlignedSegment's fields in-place
    if the result is good.
    Ideally, these records are accessed through a name-sorted AlignmentFile,
    as it throws an index error.
    Might fix this later, but I don't like it very much. If there were cython
    bindings to bwa, that would be the perfect use of it. Ehhh..
    """
    cdef int lAlignedArr, lbf
    cdef float af
    cdef cystr FastqStr, commandStr
    cdef list alignedArr, letters, numbers, tags, tag
    if(read.opt("AF") == 1.):
        # Nothing wrong a fully-aligned read.
        return read
    try:
        assert ref != "default"
    except AssertionError:
        raise Tim("Reference required for AlnRealignAS!")
    FastqStr = ASToFastqSingle(read)
    commandStr = "echo \'%s\' | bwa bwasw -M %s -" % (FastqStr, ref)
    printlog("Command String: %s" % commandStr)
    try:
        alignedArr = [i.split("\t") for i in
                      check_output(commandStr, shell=True).split("\n")
                      if not(i == "" or i[0] == "@")]
    except subprocess.CalledProcessError:
        pl("Problem with read - usually the quality scores have a "
           "quotation mark therein, which can make a realignment infeasible."
           " Returning original read.", level=logging.INFO)
        return read
    alignedArr = [i for i in alignedArr if
                  not (int(i[1]) & 256 or int(i[1]) & 2048)]
    lAlignedArr = len(alignedArr)
    if(lAlignedArr == 0):
        printlog("Could not find a suitable alignment with bwasw",
                 level=logging.INFO)
        return read
    if(lAlignedArr > 1):
        printlog("Somehow there are multiple primary alignments. I don't "
                 "know what to do, so just returning the original "
                 "read.", level=logging.INFO)
        return read
    else:
        alignedArr = alignedArr[0]
    if(alignedArr[4] == "0"):
        printlog("bwasw could not assign a non-0 MQ alignment. "
                 "Returning original read.", level=logging.INFO)
        return read
    # Now check the cigar string.
    ra = regex_compile("[A-Z]")
    rn = regex_compile("[0-9]")
    letters = [i for i in rn.split(alignedArr[5]) if i != ""]
    numbers = [int(i) for i in ra.split(alignedArr[5])[:-1]]
    # Calculate aligned fraction. Fail realignment if minAF not met.
    af = sum([i[1] for i in zip(letters, numbers)
              if i[0] == "M"]) / (1. * sum(numbers))
    if(af < minAF):
        printlog("bwasw could not align more than %s. (aligned: " % minAF +
                 "%s) Returning original read." % af, level=logging.INFO)
        return read
    # Get the realignment tags and overwrite original tags with them.
    tags = [i.split(":")[::2] for i in alignedArr[11:]]  # Gets tags
    for tag in tags:
        try:
            tag[1] = int(tag[1])
            read.set_tag(tag[0], tag[1], "i")
            continue
        except ValueError:
            pass
        try:
            tag[1] = float(tag[1])
            read.set_tag(tag[0], tag[1], "f")
            continue
        except ValueError:
            pass
        read.set_tag(tag[0], tag[1])
    read.cigarstring = alignedArr[5]
    read.reference_id = alignmentfileObj.gettid(alignedArr[2])
    read.pos = int(alignedArr[3])
    read.mapq = int(alignedArr[4])
    read.set_tag("RA", "bwasw", "Z")
    return read


@cython.returns(dict)
def makeinfodict(pysam.ctabixproxies.VCFProxy rec):
    """
    Returns a dictionary of info fields for a tabix VCF Proxy
    """
    return dict([i.split("=") for i in rec.info.split(";")])


@cython.returns(dict)
def makeformatdict(pysam.ctabixproxies.VCFProxy rec):
    """
    Returns a dictionary of format fields for a tabix VCF Proxy
    """
    return dict(zip(rec.format.split(":"), rec[0].split(":")))


@cython.returns(bint)
def DeaminationConfTest(pysam.ctabixproxies.VCFProxy rec,
                        np.longdouble_t ctfreq=-1., np.longdouble_t conf=1e-3):
    """
    Tests whether or not a given record should pass.
    """
    cdef dict InfoDict
    cdef dict Counts
    cdef int ACR
    cdef int ceiling
    cdef int AC
    InfoDict = makeinfodict(rec)
    Counts = dict([i.split(">") for i in
                  InfoDict["MACS"].split(",")])
    cons = InfoDict["CONS"]
    if(cons == "T" or cons == "A" or rec.alt == "C" or rec.alt == "G"):
        return True
    if(rec.alt == "T"):
        if(rec.ref != "C" and InfoDict["CONS"] != "C"):
            return True
        ACR = int(Counts["C"])
        AC = int(Counts["T"])
    elif(rec.alt == "A"):
        if(rec.ref != "G" and InfoDict["CONS"] != "G"):
            return True
        ACR = int(Counts["G"])
        AC = int(Counts["A"])
    if(ACR == 0):
        return True

    try:
        assert ctfreq >= 0
    except AssertionError:
        raise Tim("kwarg ctfreq required for DeaminationConfTest!")
    ceiling = GetCeiling(ACR, p=ctfreq, pVal=conf)
    if(AC > ceiling):
        return True
    return False


def PartialDeaminationConfTest(np.longdouble_t ctfreq,
                               np.longdouble_t conf=1e-3):
    """
    Returns a function that can be easily used by AbstractVCFProxyFilter.
    """
    return partial(DeaminationConfTest, ctfreq=ctfreq, conf=conf)


class AbstractVCFProxyFilter(object):
    """
    Abstract class which serves as a template for VCF record post-filtering.
    """
    def __init__(self, filterStr, func=Tim,
                 key="default", value="*"):
        if(func == Tim):
            func("func must be overridden for "
                 "AbstractVCFProxyFilter to work!")
        if(filterStr == "default"):
            raise Tim("filterStr must be set to append or replace the "
                      "filter field.")
        if(hasattr(func, "__call__")):
            self.func = func
        else:
            raise AttributeError("func must be callable!")
        self.filterStr = filterStr
        self.key = key
        self.value = value

    @cython.returns(pysam.ctabixproxies.VCFProxy)
    def __call__(self, pysam.ctabixproxies.VCFProxy rec):
        """
        You can specify multiple info tags with their values
        so long as you have comma-separated fields of equal length
        between "key" and "value".
        """
        if(not self.func(rec)):
            if(rec.filter == "PASS"):
                rec.filter = self.filterStr
            else:
                rec.filter += ";" + self.filterStr
        if(self.key != "default"):
            if("," not in self.key):
                rec.info += ";" + self.key + "=" + self.value
            else:
                for subkey, subvalue in zip(
                        self.key.split(","),
                        self.value.split(",")):
                    rec.info += ";" + subkey + "=" + subvalue
        return rec


def MakeVCFProxyDeaminationFilter(np.longdouble_t ctfreq,
                                  np.longdouble_t conf=1e-3,
                                  key="default", value="*"):
    """
    Returns the VCFProxyFilter object I wanted.
    """
    return AbstractVCFProxyFilter("DeaminationNoise",
                                  func=PartialDeaminationConfTest(ctfreq,
                                                                  conf=conf),
                                  key=key, value=value)


def SortBgzipAndTabixVCF(inVCF, outVCF="default", vcflib=True):
    """
    Sorts and tabix indexes a VCF.
    If vcflib is True, use vcfstreamsort from ekg's excellent vcf API.
    If not, hack it together with simple shell calls.
    """
    if(outVCF == "default"):
        outVCF = TrimExt(inVCF) + ".sort.vcf"
    if(not vcflib):
        if(outVCF[-3:] == ".gz"):
            check_call("zcat %s | head -n 2000 | grep '^#' > %s" % (inVCF,
                                                                    outVCF),
                       shell=True)
            check_call("zcat %s | grep -v '^#' | sort " % inVCF +
                       "-k1,1 -k2,2n >> %s" % outVCF, shell=True)
        else:
            check_call("cat %s | head -n 2000 | grep '^#' > %s" % (inVCF,
                                                                   outVCF),
                       shell=True)
            check_call("cat %s | grep -v '^#' | sort " % inVCF +
                       "-k1,1 -k2,2n >> %s" % outVCF,
                       shell=True)
    else:
        if(inVCF[-3:] == ".gz"):
            check_call("zcat %s | vcfstreamsort -w 100000 - > %s" % (inVCF,
                                                                     outVCF),
                       shell=True)
        else:
            check_call("vcfstreamsort -w 100000 %s > %s" % (inVCF, outVCF),
                       shell=True)
    check_call(["bgzip", outVCF])
    check_call(["tabix", outVCF + ".gz"])
    return outVCF + ".gz"


@cython.returns(cystr)
def SplitBed(cystr bedpath):
    bedbase = ".".join(bedpath.split(".")[0:-1])
    bedlines = ParseBed(bedpath)
    contigs = list(set(map(oig0, bedlines)))
    contigListSets = [[line for line in bedlines if line[0] == contig]
                      for contig in contigs]
    for contig, lset in zip(contigs, contigListSets):
        for line in lset:
            line[1] -= 20
            line[2] += 20
        open(bedbase + ".split." + contig + ".bed", "w").writelines(
            "\n".join(["\t".join(map(str, arr)) for arr in lset]))
    return ":".join([bedbase + ".split." + contig + ".bed" for
                     contig in contigs])


@cython.returns(cystr)
def MergeBamList(bamlist, picardpath="default", memStr="-Xmx6G",
                 outbam="default"):
    """
    Merges a list of BAMs. Used for merging discordant read bams for
    parallelized VCF calls.
    """
    if(isinstance(bamlist, str)):
        bamlist = bamlist.split(":")
    if(outbam == "default"):
        outbam = bamlist[0].split(".")[0:-1] + ".merged.bam"
    commandStr = "java %s -jar %s AS=true" % (memStr, picardpath)
    commandStr += " I=" + " I=".join(bamlist)
    commandStr += " O=%s" % outbam
    check_call(shlex.split(commandStr))
    return outbam


def DevNullPopen(string):
    return subprocess.Popen(string, shell=True, stdout=open(os.devnull, 'w'))


class PopenCall(object):
    """
    Contains a Popen call and the command string used.
    """
    def __init__(self, string, maxresubs=10):
        self.commandString = string
        self.popen = DevNullPopen(string)
        self.poll = self.popen.poll
        self.communicate = self.popen.communicate
        self.resubmissions = 0
        self.maxresubs = maxresubs

    def resubmit(self):
        if(self.resubmissions >= self.maxresubs):
            raise Tim(
                "Resubmission limit %s reached!" % self.maxresubs)
        self.popen = DevNullPopen(self.commandString)
        self.resubmissions += 1


@cython.returns(cystr)
def GetOutVCFFromBMFsnvPopen(d):
    return d.commandString.split(" ")[9]


@cython.returns(cystr)
def GetOutVCFFromBMFsnvCStr(d):
    return d.split(" ")[9]


def GetOutBAMForSamtoolsViewCStr(d):
    return d.split(";")[0].split(" ")[-2]


def ReturnDefault(x):
    return "default"


class PopenDispatcher(object):
    """
    Contains a list of strings for Popen to use and a set of Popen options,
    and a list of threads it should have going at once.
    Cleanup must be a function to call at the end of daemon's execution.
    """
    defaultSleeptime = 5

    def __init__(self, stringlist, threads=4, MaxResubmissions=50,
                 func=None, cleanup=None, sleeptime=defaultSleeptime):
        print("Initializing PopenDispatcher!")
        assert len(stringlist) > 0
        self.queue = deque(stringlist)
        self.dispatches = []
        self.threadcount = threads - 1  # Main thread counts as one, too.
        self.outstrs = {cStr: None for cStr in stringlist}
        self.completed = 0
        self.jobnumber = 1
        self.submitted = 0
        self.alljobssubmitted = False
        self.alljobscompleted = False
        self.resubmittedjobcounts = 0
        self.MaxResubmissions = MaxResubmissions
        if(hasattr(func, "__call__") is False):
            raise Tim("func must be callable!")
        self.getReturnValue = func
        self.bgStrs = []
        self.fgStrs = []
        if(cleanup is not None):
            self.cleanup = cleanup
        self.sleeptime = sleeptime

    @cython.returns(int)
    def _getJobNumber(self):
        return self.submitted + 1

    def _submit(self):
        """
        Function for submitting bg (background) jobs.
        """
        if(len(self.queue) == 0):
            print("All jobs already submitted. :)")
            return
        while(len(self.dispatches) < self.threadcount):
            if(len(self.queue) != 0):
                print("Submitting job number %s." % str(self._getJobNumber()))
                self.submitted += 1
                cStr = self.queue.popleft()
                self.bgStrs.append(cStr)
                submitted = PopenCall(cStr)
                print("Command String: %s" % submitted.commandString)
                self.dispatches.append(submitted)
            else:
                #print("All jobs submitted - check in later.")
                pass

    def _check(self):
        """
        Checks for finished processes. If they are, knock them off the list.
        """
        for d in self.dispatches:
            if(d.poll() is not None):
                if(d.poll() != 0):
                    # See if calling it again solves the problem
                    print("Had a non-zero return status. Trying again!")
                    d.resubmit()
                    self.resubmittedjobcounts += 1
                    if(self.resubmittedjobcounts > self.MaxResubmissions):
                        raise CalledProcessError(
                            d.poll(), d.commandString,
                            "So we resubmitted jobs until we passed the "
                            "limit for failed jobs. "
                            "(%s). Oops!" % self.MaxResubmissions)
                    continue
                self.completed += 1
                self.outstrs[
                    d.commandString] = self.getReturnValue(d.commandString)
                self.dispatches.remove(d)
        return len(self.dispatches)

    def daemon(self):
        """
        Little daemon hides beyond sight of the world, looking for mischief
        and maintaining order. We don't have to watch, but he's always
        lurking, always working, in the night never shirking.
        He hasn't slept since the dawn of time - unweary, inhuman, lost in the
        world that was wild and waste.

        It also submits foreground jobs.
        """
        self._submit()
        while len(self.queue) != 0:
            print("Submitting set of jobs for daemon.")
            self._submit()
            newqueue = tee(self.queue, 1)[0]
            try:
                nextStr = newqueue.next()
            except StopIteration:
                time.sleep(self.sleeptime)
                continue
            if("#" in nextStr):
                fgCStr = self.queue.popleft()
                self.fgStrs.append(fgCStr)
                print("Foreground submitting job #%s" % self._getJobNumber())
                try:
                    exec(fgCStr.split("#")[1])
                    self.submitted += 1
                except Exception:
                    print("Failed function call: %s" % fgCStr.split("#")[1])
                    raise FunctionCallException(
                        fgCStr, "Foreground eval call failed.",
                        False)
                self.outstrs[fgCStr] = self.getReturnValue(fgCStr)
            else:
                time.sleep(self.sleeptime)
            self._check()
        print("All jobs submitted! Yay.")
        while(self._check() > 0 and len(self.queue) != 0):
            time.sleep(self.sleeptime)
            threadcount = self._check()
            if(threadcount < self.threadcount and len(self.queue) != 0):
                self._submit()
        while(self._check() != 0):
            """
            while(sum([d.popen.poll() is not None for d in
                       self.dispatches]) < len(self.dispatches)):
                pl("Sleeping because some jobs haven't returned yet.")
                time.sleep(self.sleeptime)
            """
            pl("Sleeping because some jobs haven't returned yet.")
            time.sleep(self.sleeptime)
        for key in self.outstrs.iterkeys():
            if(self.outstrs[key] is None):
                print("fgStrs: %s" % ":".join(self.fgStrs))
                print("bgStrs: %s" % ":".join(self.bgStrs))
                print(repr(self.outstrs))
                raise Tim("Command string %s has a key value of None!" % key)
        if("cleanup" in dir(self)):
            self.cleanup()
        return 0


def GetBamBedList(bampath, bedpath):
    bedlist = SplitBed(bedpath).split(":")
    bamlist = [".".join(bampath.split(".")[0:-1]) + "." +
               i.split(".")[-2] + ".bam" for i in bedlist]
    return zip(bedlist, bamlist)


def GetSplitBAMCStrs(bampath, ziplist):
    return ["samtools view -bh -L "
            "%s -o %s %s;samtools index %s" % (bed, bam, bampath, bam)
            for bed, bam in ziplist]


def GetSplitBamPopen(bampath, bedpath, threads=4):
    ziplist = GetBamBedList(bampath, bedpath)
    cStrs = GetSplitBAMCStrs(bampath, ziplist)
    return PopenDispatcher(cStrs, threads=threads,
                           func=GetOutBAMForSamtoolsViewCStr)


def SplitBamParallel(bampath, bedpath, threads=4):
    """
    Creates a PopenDispatcher object and calls them all in parallel.
    """
    pl("Splitting BAM %s by BED %s" % (bampath, bedpath))
    PD = GetSplitBamPopen(bampath, bedpath, threads=threads)
    return PD.daemon()


def SplitBamByBed(bampath, bedpath):
    ziplist = GetBamBedList(bampath, bedpath)
    for bed, bam in ziplist:
        # Get the reads that overlap region of interest.
        check_call(["samtools", "view", "-bh", "-L", bed,
                    "-o", bam, bampath])
        # Index it so that variant-calling can happen.
        check_call(["samtools", "index", bam])
        if(not ospath.isfile(bam + ".bai")):
            raise Tim("Index for bam not created!")
    return ziplist


def SplitBamByBedPysam(bampath, bedpath):
    """
    Rather than standard split by bed which uses parallel shell calls,
    instead this does it manually via pysam. Not sure if it's faster or not.
    """
    cdef AlignedSegment_t rec
    cdef pysam.calignmentfile.AlignmentFile inHandle
    cdef cystr bam
    pl("Getting bamlist.")
    bamlist = map(oig1, GetBamBedList(bampath, bedpath))
    refContigNumList = [PysamToChrDict[bam.split(".")[-2]] for bam in bamlist]
    inHandle = pysam.AlignmentFile(bampath, "rb")
    handles = [pysam.AlignmentFile(bam, "wb", template=inHandle) for
               bam in bamlist]
    refHandleMap = dict(zip(refContigNumList, handles))
    # print(repr(refHandleMap))
    pl("Now splitting bam into one bam per contig.")
    for rec in inHandle:
        if(rec.is_unmapped):
            continue
        try:
            refHandleMap[rec.reference_id].write(rec)
        except KeyError:
            continue
    for bam in bamlist:
        check_call(["samtools", "index", bam])
        if(ospath.isfile(bam + ".bai") is False):
            raise Tim("samtools couldn't index this bam. "
                      "It should already be sorted. Abort!")
    return bamlist


def SlaveDMP(bsFastq1, bsFastq2,
             p3Seq="default", p5Seq="default",
             overlapLen=6, sortMem=None, head=None):
    from MawCluster import BCFastq
    outfq1 = TrimExt(bsFastq1) + ".dmp.fastq"
    outfq2 = TrimExt(bsFastq2) + ".dmp.fastq"
    if(sortMem is None):
        sortMem = "768M"
    if(head is None):
        head = 4
    sortFastq1, sortFastq2 = BCFastq.BarcodeSortBoth(bsFastq1, bsFastq2,
                                                     sortMem=sortMem)
    consFastq1, consFastq2 = BCFastq.pairedFastqConsolidate(
        sortFastq1, sortFastq2)
    trimFastq1, trimFastq2 = BCFastq.CutadaptPaired(
            consFastq1, consFastq2, overlapLen=overlapLen,
            p3Seq=p3Seq, p5Seq=p5Seq, outfq1=outfq1, outfq2=outfq2)
    return trimFastq1, trimFastq2


@cython.returns(cystr)
def SlaveDMPCommandString(cystr bsFastq1, cystr bsFastq2,
                          cystr sortMem=None,
                          overlapLen=None, head=None,
                          p3Seq=None, p5Seq=None):
    """
    Returns a command string for calling bmftools snv
    """
    if(overlapLen is None):
        overlapLen = 6
    if(p3Seq is None or p5Seq is None):
        raise UnsetRequiredParameter(
            "p3Seq and p5Seq must both be set to run SlaveDMPCommandString "
            "because cutadapt needs this information.")
    cStr = ("python -c 'from utilBMF.HTSUtils import SlaveDMP;SlaveDMP"
            "(\"%s\",\"%s\", sortMem=\"%s\"" % (bsFastq1, bsFastq2, sortMem) +
            ", overlapLen=%s, head=%s, p3Seq=\"" % (overlapLen, head) +
            "\%s\", p5Seq=\"%s\")'" % (p3Seq, p5Seq))
    FnCall = ("from utilBMF.HTSUtils import SlaveDMP;SlaveDMP("
              "\"%s\",\"%s\", sortMem=\"%s" % (bsFastq1, bsFastq2, sortMem) +
              "\", overlapLen=%s, head=%s, p3Seq=" % (overlapLen, head) +
              "\"%s\", p5Seq=\"%s\")" % (p3Seq, p5Seq))
    return cStr + " #" + FnCall


@cython.returns(cystr)
def GetFastqPathsFromDMPCStr(cystr cStr):
    """
    Helper function for a PopenDispacher for SlaveDMPCommandString.
    """
    return ",".join([i.replace("'", "").replace(
        "\"", "").replace(".fastq", ".dmp.fastq") for
                     i in cStr.split(";")[1].split("(")[1].split(",")[:2]])


def GetParallelDMPPopen(fqPairList, sortMem=None, int threads=-1,
                        head=None, overlapLen=None, p3Seq=None,
                        p5Seq=None):
    """
    Makes a PopenDispatcher object for calling these variant callers.
    """
    if(threads < 0):
        raise UnsetRequiredParameter("threads must be set for"
                                     " GetParallelDMPopen.")
    pl("Dispatching BMF dmp jobs")
    return PopenDispatcher([SlaveDMPCommandString(*fqPair, head=head,
                                                  sortMem=sortMem,
                                                  overlapLen=overlapLen,
                                                  p3Seq=p3Seq, p5Seq=p5Seq) for
                            fqPair in fqPairList],
                           threads=threads,
                           func=GetFastqPathsFromDMPCStr)


@cython.returns(cystr)
def BMFsnvCommandString(cystr bampath, cystr conf="default",
                        cystr bed="default"):
    """
    Returns a command string for calling bmftools snv
    """
    if(conf == "default"):
        raise Tim("conf must be set for BMFsnvCommandString.")
    tag = str(uuid.uuid4().get_hex().upper()[0:8])
    if(bed == "default"):
        raise Tim("bed must be set for BMFsnvCommandString.")
    outVCF = ".".join(bampath.split(".")[:-1] +
                      ["bmf", "vcf"])
    config = parseConfig(conf)
    minMQ, minBQ, minFA = config["minMQ"], config["minBQ"], config["minFA"]
    minFracAgreed, MaxPValue = config["minFracAgreed"], config["MaxPValue"]
    ref = config["ref"]
    FnCall = ("from MawCluster.VCFWriters import SNVCrawler;"
              "SNVCrawler('%s', bed='%s', minMQ=" % (bampath, bed) +
              "%s, minBQ=%s, minFA=%s, " % (minMQ, minBQ, minFA) +
              "minFracAgreed=%s, MaxPValue=" % (minFracAgreed) +
              "%s, reference='%s', OutVCF='%s'," % (MaxPValue, ref, outVCF) +
              "writeHeader=False)")
    return ("bmftools snv --conf %s %s " % (conf, bampath) +
            "--bed %s --is-slave --outVCF %s #%s" % (bed, outVCF, FnCall))


@cython.returns(cystr)
def GetUUIDFromCommandString(cystr cStr):
    """
    Gets the UUID tag from slave jobs. Used for concatenating VCFs from
    parallel calls.
    """
    return cStr.split("|")[1]


def GetBMFsnvPopen(bampath, bedpath, conf="default", threads=4,
                   parallel=False):
    """
    Makes a PopenDispatcher object for calling these variant callers.
    """
    if conf == "default":
        raise Tim("conf file must be set for GetBMFsnvPopen")
    ziplist = GetBamBedList(bampath, bedpath)
    if(parallel is True):
        SplitBamParallel(bampath, bedpath)
    else:
        SplitBamByBedPysam(bampath, bedpath)
    pl("Dispatching BMF jobs")
    return PopenDispatcher([BMFsnvCommandString(tup[1],
                                                conf=conf,
                                                bed=tup[0]) for
                            tup in ziplist],
                           threads=threads,
                           func=GetOutVCFFromBMFsnvCStr)


@cython.returns(cystr)
def TrimExt(cystr fname, cystr exclude=""):
    """
    Trims the extension from the filename so that I don't have to type this
    every single time.
    :param fname - filename to strip
    :param exclude - string for exclusion.
    (Any number of comma-separated words may be provided.)
    """
    cdef cystr tmpStr
    if(fname is not None):
        tmpList = fname.split("/")[-1].split(".")[:-1]
        if(tmpList[-1] == "gz"):
            tmpList = tmpList[:-1]
        try:
            return ".".join([tmpStr for tmpStr in tmpList if
                             tmpStr not in exclude.split(",")])
        except IndexError:
            return tmpList[0]
    else:
        raise Tim("Cannot trim an extension of a None value!")


@cython.returns(cystr)
def NPadSequence(cystr seq, int n=300):
    """
    Pads a sequence with "n" Ns.
    """
    return "N" * n + seq + "N" * n


@cython.returns(cystr)
def FastaStyleSequence(cystr seq):
    return ">" + seq + "\n" + seq


@cython.returns(cystr)
def PadAndMakeFasta(cystr seq, int n=300):
    return FastaStyleSequence(NPadSequence(seq, n=n))


@cython.returns(list)
def GetUniqueItemsD(dict inDict):
    return [i[0] for i in inDict.iteritems() if i[1] == 1]


@cython.returns(list)
def GetUniqueItemsL(list inList):
    return GetUniqueItemsD(cyfreq(inList))


def GetDSIndels(inBAM, outBAM):
    """
    Gets read pairs tagged with DS[ID], sorts by coordinate (to facilitate
    piling up the deletions)

    Note: It is faster to do samtools sort -O... then grep, then samtools view
    than to do samtools view -Sbh | grep, then samtools sort.
    And, of course, it's faster to get a conversion while sorting.
    """
    randPref = str(uuid.uuid4().get_hex().upper()[0:8]) + ".hey.i.am.a.prefix"
    cStr = ("samtools sort -O sam -T %s %s| grep 'DS[ID]\|^@' | " % (randPref,
                                                                     inBAM) +
            "samtools view -Sbh - > %s" % outBAM)
    pl("About to get DSIndels with command string: %s" % cStr)
    check_call(cStr)
    return outBAM


def hamming_cousins_exact(cystr s, int n,
                          list alphabet=["A", "C", "G", "T"]):
    """Generate strings over alphabet whose Hamming distance from s is
    exactly n.

    >>> sorted(hamming_cousins_exact('abc', 0, 'abc'))
    ['abc']
    >>> sorted(hamming_cousins_exact('abc', 1, 'abc'))
    ['aac', 'aba', 'abb', 'acc', 'bbc', 'cbc']
    >>> sorted(hamming_cousins_exact('aaa', 2, 'ab'))
    ['abb', 'bab', 'bba']

    """
    cdef int lalphabetSub1
    lalphabetSub1 = len(alphabet) - 1
    for positions in combinations(xrange(len(s)), n):
        for replacements in product(range(lalphabetSub1), repeat=n):
            cousin = list(s)
            for p, r in zip(positions, replacements):
                if cousin[p] == alphabet[r]:
                    cousin[p] = alphabet[-1]
                else:
                    cousin[p] = alphabet[r]
            yield ''.join(cousin)


def hamming_cousins(cystr s, int n=0,
                    list alphabet=["A", "C", "G", "T"]):
    """Generate strings over alphabet whose Hamming distance from s is
    less than or equal to n.

    >>> sorted(hamming_cousins('abc', 0, 'abc'))
    ['abc']
    >>> sorted(hamming_cousins('abc', 1, 'abc'))
    ['aac', 'aba', 'abb', 'abc', 'acc', 'bbc', 'cbc']
    >>> sorted(hamming_cousins('aaa', 2, 'ab'))
    ['aaa', 'aab', 'aba', 'abb', 'baa', 'bab', 'bba']

    """
    cdef int i
    return chain(*(hamming_cousins_exact(s, i, alphabet) for i
                   in range(n + 1)))


@cython.returns(Insertion_t)
def GetInsertionFromAlignedSegment(AlignedSegment_t read,
                                   pysam.cfaidx.FastaFile handle=None):
    """
    Creates an Insertion object under the assumption that there is
    only one set of inserted bases in the read. I should come back and
    generalize it...
    """
    cdef tuple InsInf
    InsInf = GetInsertedStrs(read)[0]
    return Insertion(read, contig=PysamToChrDict[read.reference_id],
                     start=InsInf[1], seq=InsInf[0], handle=handle)


@cython.returns(Deletion_t)
def GetDeletionFromAlignedSegment(AlignedSegment_t read,
                                  pysam.cfaidx.FastaFile handle=None):
    """
    Creates a Deletion object from a read.
    """
    cdef int start, end
    cdef list coords
    coords = GetDeletedCoordinates(read)
    start = coords[0]
    end = coords[-1]
    return Deletion(read, contig=PysamToChrDict[read.reference_id],
                    start=start, end=end, handle=handle)


@cython.returns(cystr)
def is_reverse_to_str(bint boolean):
    if(boolean):
        return "reverse"
    elif(boolean is False):
        return "forward"
    else:
        return "unmapped"


@cython.returns(cystr)
def ssStringFromRead(AlignedSegment_t read):
    return ("#".join(map(str, sorted([read.reference_start,
                                      read.reference_end]))) +
            "#%s" % (is_reverse_to_str(read)))


cdef class AbstractIndelContainer(object):
    """
    Base class for insertion and deletion container objects.

    Type can be -1, 0, or 1. (Deletion, deletion and insertion, and
    just insertion)
    Start and end refer to different things for insertions and deletions.
    For a deletion, start is the first missing base and end is the last missing
    reference base position.
    seq should be None for a deletion
    """

    def __init__(self, cystr contig, int start=-666,
                 int end=-1, int type=-137,
                 cystr seq=None):
        self.contig = contig
        self.start = start
        self.end = end
        self.type = type
        self.seq = seq
        self.readnames = []
        self.uniqStr = None
        self.StartStops = []

    def __str__(self):
        raise Tim("Abstract method must be inherited. Sorry, cdef won't let "
                  "me actually make this an abstract class.")

    @cython.returns(int)
    def __len__(self):
        """
        Returns the number of reads supporting it which have been queried
        against this object since its creation.
        """
        return len(self.readnames)

    @cython.returns(bint)
    def compare(self, indelObj):
        """
        Used for comparing two insertion or deletion objects.
        This just makes it a little cleaner and more flexible.
        That way, if they have the same unique identifier string
        but a different number of readnames.
        """
        assert isinstance(indelObj, AbstractIndelContainer)
        return self.uniqStr == indelObj.uniqStr

    def add(self, indelObj, inplace=True):
        """
        If the uniqStr attributes are identical, merge the families.
        I imagine that this would be useful when comparing lots of sets
        of indels across different samples or doing multiple passes.
        (e.g., duplex-supported indels and otherwise)
        """
        assert isinstance(indelObj, AbstractIndelContainer)
        if self.compare(indelObj):
            if(inplace):
                self.readnames += indelObj.readnames
                return self
            else:
                newIndelObj = ccopy(indelObj)
                newIndelObj.readnames += self.readnames
                return newIndelObj

    @cython.returns(dict)
    def GetNameCounter(self):
        return cyfreq(self.readnames)

    @cython.returns(cystr)
    def __getitem__(self, int index):
        return self.readnames[index]

    @cython.returns(list)
    def sort(self):
        self.readnames = sorted(self.readnames)

    @cython.returns(int)
    def getNumSS(self):
        return(len(set(self.StartStops)))

    def register(self, AlignedSegment_t read):
        self.readnames.append(read.query_name)
        self.StartStops.append(ssStringFromRead(read))

    def merge(self, AbstractIndelContainer_t AIC):
        try:
            assert AIC.uniqStr == self.uniqStr
        except AssertionError:
            raise Tim("To merge two IndelContainer objects, "
                      "their unique string description must match.")
        self.readnames += AIC.readnames
        self.StartStops += AIC.StartStops


cdef class Insertion(AbstractIndelContainer):
    """
    Expansion of the AbstractIndelContainer object.
    "Start" is actually the preceding base to the insertion (counting up),
    "End" is the following. end should always be greater than start.
    The important thing is to be able to hold read names and have a unique
    string representing each indel so that we can make calls.

    If no handle is provided, shen (Shannon Entropy) is set to -1.
    """

    def __init__(self, AlignedSegment_t read,
                 cystr contig, int start=-1,
                 cystr seq=None, pysam.cfaidx.FastaFile handle=None,
                 int window=20):
        if(start < 0):
            raise Tim("start required for InsertionContainer.")
        self.contig = contig
        self.start = start
        self.end = start + 1
        self.type = 1
        self.seq = seq
        self.readnames = [read.query_name]
        self.uniqStr = "Insertion|%s:%s,%s|%s" % (self.contig, self.start,
                                                  self.end, self.seq)
        try:
            self.shen = min([shen(handle.fetch(read.reference_id,
                                               start=start - window,
                                               end=start)),
                             shen(handle.fetch(read.reference_id,
                                               start=self.end,
                                               end=self.end + window))])
        except AttributeError:
            self.shen = -1.
        self.handle = handle
        self.shenwindow = window
        self.StartStops = [ssStringFromRead(read)]

    @cython.returns(cystr)
    def __str__(self):
        return self.uniqStr + "|%s" % len(self.readnames)


cdef class Deletion(AbstractIndelContainer):
    """
    Expansion of the AbstractIndelContainer object.
    "Start" is the first missing base.
    "End" is the last missing base.
    If the deletion is of length one, start and end should be the same.
    """

    def __init__(self, AlignedSegment_t read,
                 cystr contig=None, int start=-1,
                 int end=-1,
                 pysam.cfaidx.FastaFile handle=None, int window=20):
        self.contig = contig
        self.start = start
        self.end = end
        self.type = -1
        self.readnames = [read.query_name]
        self.uniqStr = "Deletion|%s:%s,%s" % (self.contig, self.start,
                                              self.end)
        try:
            self.shen = min([shen(handle.fetch(contig,
                                               start=start - window,
                                               end=start)),
                             shen(handle.fetch(contig,
                                               start=end, end=end + window))])
        except AttributeError:
            self.shen = -1.
        self.handle = handle
        self.shenwindow = window
        self.StartStops = [ssStringFromRead(read)]

    @cython.returns(cystr)
    def __str__(self):
        return self.uniqStr + "|%s" % len(self.readnames)


cdef class IndelQuiver(object):

    """
    Class for holding on to the indel objects. Holds a list for each
    type of indel, as well as a dictionary where the keys are unique
    string descriptors for the mutation and the values are lists
    of read names for reads supporting that mutation.
    Counts is a similar object, but with the length of the data
    field as a value instead of the list itself.
    """
    def __init__(self, cystr ref=None, int window=10,
                 int minMQ=0, int minFM=0,
                 cystr bam=None, int minNumSS=0,
                 float minShen=0.2, int minPairs=1):
        self.data = {}
        self.readnames = {}
        self.counts = {}
        self.fastaRef = pysam.FastaFile(ref)
        self.bam = pysam.AlignmentFile(bam, "rb")
        self.window = window
        self.minMQ = minMQ
        self.minFM = minFM
        self.minPairs = minPairs
        self.minShen = minShen
        self.minNumSS = minNumSS

    @cython.returns(int)
    def __len__(self):
        return len(self.data)

    @cython.returns(list)
    def __getitem__(self, cystr key):
        return self.data[key]

    def __setitem__(self, cystr key, list value):
        self.data[key] = value

    def setIndelShen(self, AbstractIndelContainer_t indelObj):
        indelObj.shen = min([shen(self.fastaRef(indelObj.contig,
                                                indelObj.start - self.window,
                                                indelObj.start)),
                             shen(self.fastaRef(indelObj.contig, indelObj.end,
                                                indelObj.end + self.window))])
        indelObj.shenwindow = self.window

    def iterkeys(self):
        return self.data.iterkeys()

    @cython.returns(list)
    def keys(self):
        return self.data.keys()

    def iteritems(self):
        return self.data.iteritems()

    @cython.returns(list)
    def items(self):
        return self.data.items()

    def itervalues(self):
        return self.data.itervalues()

    @cython.returns(list)
    def values(self):
        return self.data.values()

    def addRead(self, AlignedSegment_t read):
        cdef cystr SVTag
        cdef AbstractIndelContainer_t Indel
        if(read.opt("FM") < self.minFM):
            return
        if(read.mapping_quality < self.minMQ):
            return
        try:
            SVTag = read.opt("SV")
        except KeyError:
            raise Tim("read must have an SV tag for IndelQuiver!")
        if("DSI" in SVTag):
            Indel = GetInsertionFromAlignedSegment(read, handle=self.fastaRef)
            try:
                self[Indel.uniqStr].register(read)
            except KeyError:
                self.setIndelShen(Indel)
                self[Indel.uniqStr] = Indel
        if("DSD" in SVTag):
            Indel = GetDeletionFromAlignedSegment(read, handle=self.fastaRef)
            try:
                self[Indel.uniqStr].register(read)
            except KeyError:
                self.setIndelShen(Indel)
                self[Indel.uniqStr] = Indel
        self.counts = {key: len(values) for key, values in
                       self.readnames.itervalues()}

    def mergeQuiver(self, IndelQuiver_t quiverObj):
        cdef cystr key
        for key in quiverObj.readnames.keys():
            try:
                self.readnames[key] += quiverObj[key]
            except KeyError:
                self.readnames[key] = quiverObj[key]
        self.counts = {key: len(values) for key, values in
                       self.readnames.itervalues()}
        for key in quiverObj.keys():
            try:
                self[key].merge(quiverObj[key])
            except KeyError:
                self[key] = quiverObj[key]

    @cython.returns(dict)
    def getIndelCounter(self, cystr uniqStr):
        try:
            return cyfreq(self.readnames[uniqStr])
        except KeyError:
            return {}

    @cython.returns(IDVCFLine_t)
    def makeVCFLine(self, AbstractIndelContainer_t IC):
        return IDVCFLine(IC, self)


cdef class BamTag(object):
    """
    Contains a tag, a value, and a type, all of which are string objects.
    """
    @classmethod
    def fromstring(cls, cystr tagString):
        """The dictionary is a pythonic proxy for a switch statement.
        Initializes and returns a BamTag object.
        """
        cdef list tokens
        tokens = tagString.split(":")
        return cls(tokens[0], tokens[1],
                   {"Z": tokens[2], "A": tokens[2], "i": int(tokens[2]),
                    "f": float(tokens[2]), "H": tokens[2],
                    "B": tokens[2]}[tokens[1]])

    @classmethod
    def fromtuple(cls, tuple tag):
        cdef cystr tagtype
        try:
            tagtype = TagTypeDict[tag[0]]
        except KeyError:
            if(isinstance(tag[1], int)):
                tagtype = "i"
            elif(isinstance(tag[1], float)):
                tagtype = "f"
            else:
                tagtype = "Z"  # A safer fallback
        return cls(tag[0], tagtype, tag[1])

    def __init__(self, cystr tag, cystr tagtype, value):
        self.tag = tag
        self.tagtype = tagtype
        self.value = value

    @cython.returns(cystr)
    def __str__(self):
        """
        In [14]: %timeit c = "PV" + ":" + "Z" + ":%s" % 1337
        10000000 loops, best of 3: 55.7 ns per loop

        In [15]: %timeit c = "PV" + ":" + "Z" + ":" +  str(1337)
        1000000 loops, best of 3: 181 ns per loop

        In [17]: %timeit c = "%s:%s:%s" % ("PV", "Z", 1337)
        1000000 loops, best of 3: 270 ns per loop
        In [19]: %timeit c = ":".join(map(str, ["PV", "Z", 1337]))
        1000000 loops, best of 3: 710 ns per loop
        """
        cdef cystr ret
        if(self.tagtype == "B"):
            if(isinstance(self.value, float)):
                ret = self.tag + ":B:f%s" % ",".join(map(str, self.value))
            else:
                ret = self.tag + ":B:i%s" % ",".join(map(str, self.value))
            return ret if(len(self.value) > 1) else ret + ","

        else:
            return self.tag + ":" + self.tagtype + ":%s" % (self.value)


cdef class IDVCFLine(object):

    """
    Contains data for writing to a VCF Line, where each ALT has its own line.

    minNumSS is the minumum number of start/stop combinations required to
    support a variant call.
    """
    def __init__(self, AbstractIndelContainer_t IC, IndelQuiver_t quiver=None):
        cdef pysam.calignedsegment.PileupColumn PileupCol
        cdef pysam.calignmentfile.IteratorColumnRegion pileupIt
        cdef int tmpCov
        cdef float MDP
        cdef tuple i
        cdef cystr key
        cdef list ffkeys
        self.CHROM = IC.contig
        self.NumStartStops = len(set(IC.StartStops))
        if(isinstance(IC, Insertion)):
            """
            I would subtract one to get the start, but I would also
            add one to correct for 0 vs 1-based.
            One is subtracted from the fetch.
            """
            self.TYPE = "ins"
            self.POS = IC.start
            self.REF = quiver.fastaRef.fetch(IC.contig,
                                             self.POS - 1,
                                             self.POS)
            self.ALT = self.REF + IC.seq
            self.LEN = len(IC.seq)
        elif(isinstance(IC, Deletion)):
            """
            Because the start of a Deletion object is the first deleted
            base, the start of the variant is that minus one.
            Then, we're incrementing for 0 vs 1
            The end of the fetch is incremented by one for 0 vs 1 and
            once more to get the base after the last deleted base.
            """
            self.TYPE = "del"
            self.POS = IC.start
            self.REF = quiver.fastaRef.fetch(IC.contig,
                                             self.POS - 1,
                                             IC.end + 1)
            self.ALT = self.REF[0]
            self.LEN = IC.end - IC.start
        else:
            raise Tim("Sorry, I haven't finished this VCF writer"
                      " for complex indels.")
        self.ID = IC.uniqStr
        self.reverseStrandFraction = sum(
            ["reverse" in ssString for ssString in
             IC.StartStops]) / len(IC.StartStops)
        self.BothStrandSupport = (self.reverseStrandFraction > 0.001 and
                                  self.reverseStrandFraction < 0.999)
        self.QUAL = 20. * len(IC)  # 20 for each supporting read. Will change!
        pileupIt = quiver.bam.pileup(self.CHROM, self.POS, max_depth=200000)
        PileupCol = pileupIt.next()
        while PileupCol.pos < self.POS - 5:
            PileupCol = pileupIt.next()
        tmpCov = 0
        while PileupCol.pos < self.POS + 5:
            tmpCov += PileupCol.n
            PileupCol = pileupIt.next()
        self.MDP = tmpCov * 1. / 10
        self.FILTER = None
        self.NDPS = len([i for i in cyfreq(IC.readnames).iteritems() if
                         i[1] > 1])
        self.DPA = len(IC.readnames)
        if(len(set(IC.readnames)) < quiver.minPairs):
            try:
                self.FILTER.append(";InsufficientReadPairs")
            except AttributeError:
                self.FILTER = "InsufficientReadPairs"
        if(self.NumStartStops < quiver.minNumSS):
            try:
                self.FILTER.append(";InsufficientStopStarts")
            except:
                self.FILTER = "InsufficientStopStarts"
        if(IC.shen < quiver.minShen):
            try:
                self.FILTER.append(";LowComplexity")
            except:
                self.FILTER = "LowComplexity"
        if(self.FILTER is None):
            self.FILTER = "PASS"
        self.InfoFields = {"SHENWINDOW": quiver.window,
                           "MINSHEN": quiver.minShen}
        self.InfoStr = ";".join([key + "=" + str(self.InfoFields[key])
                                 for key in sorted(self.InfoFields.keys())])
        self.FormatFields = {"SHEN": IC.shen, "TYPE": self.TYPE,
                             "NDPS": self.NDPS, "DPA": self.DPA,
                             "MDP": self.MDP, "LEN": self.LEN}
        ffkeys = sorted(self.FormatFields.keys())
        self.FormatStr = (
            ":".join(ffkeys) +
            "\t" + ":".join(str(
                self.FormatFields[key]) for key in ffkeys))
        self.str = "\t".join(map(str, [self.CHROM,
                                       self.POS, self.ID,
                                       self.CONS, self.ALT,
                                       self.QUAL, self.FILTER,
                                       self.InfoStr, self.FormatStr]))

    def update(self):
        cdef list ffkeys
        ffkeys = sorted(self.FormatFields.keys())
        self.FormatStr = (":".join(ffkeys) + "\t" +
                          ":".join([str(self.FormatFields[key]) for
                                    key in ffkeys]))
        self.InfoStr = ";".join([key + "=" + str(self.InfoFields[key])
                                 for key in sorted(self.InfoFields.keys())])

    def __str__(self):
        self.update()
        self.str = "\t".join(map(str, [self.CHROM, self.POS,
                                       self.ID, self.REF, self.ALT,
                                       self.QUAL, self.FILTER, self.InfoStr,
                                       self.FormatStr]))
        return self.str


@cython.returns(list)
def GetSCFractionArray(inBAM):
    cdef pysam.calignmentfile.AlignmentFile inHandle
    cdef AlignedSegment_t i
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    return [FractionSoftClipped(i) for i in inHandle]


@cython.returns(ndarray)
def GetTlenArray(inBAM):
    cdef AlignedSegment_t i
    cdef ndarray[np.int64_t, ndim=2] tlens
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    tlens = np.array([i.tlen for i in inHandle if i.is_read1 and i.tlen != 0],
                     dtype=np.int64, ndmin=2)
    return np.absolute(tlens)


def PlotTlen(inBAM, outfile="default"):
    cdef ndarray[np.int64_t, ndim=2] tlens
    pl("About to load tlens")
    tlens = GetTlenArray(inBAM)
    pl("Successfully loaded tlens.")
    import matplotlib as mpl
    mpl.use('Agg')      # With this line = figure disappears
    import matplotlib.pyplot as plt
    import matplotlib.mlab as mlab
    tlens = tlens.reshape(-1, 1)
    pl("Looks like I reshaped successfully.")
    mu, sigma = nmean(tlens[tlens < 5000]), cyStdFlt(tlens[tlens < 5000])
    # Only includes them in that calculation for reasonably possible
    # insert lengths.
    n, bins, patches = plt.hist(tlens, 200, normed=1, facecolor='green',
                                alpha=0.75)
    y = mlab.normpdf(bins, mu, sigma)
    l = plt.plot(bins, y, 'r--', linewidth=1)
    plt.xlabel('Number of read pairs with inferred template length')
    plt.ylabel('Inferred template length')
    plt.title(r'$\mathrm{Histogram\ of\ inferred template lengths:}\ '
              r'\mu=%s,\ \sigma=%s$' % (mu, sigma))
    plt.axis([nmin(tlens), nmax(tlens[tlens < 5000]), 0., nmax(n)])
    plt.grid(True)
    plt.savefig(outfile + ".png")
    return outfile


def itersplit(inString, regexStr="\w+"):
    """
    Returns a split iterator over a string with a given regex.
    Default regexStr causes it to iterate over words.
    """
    return (map(mcgroup, finditer(regexStr, inString)))


@cython.returns(set)
def GetBamTagTypes(cystr bamfilestring):
    """
    :param bamfilestring: full text of a bam file
    """
    return set(itersplit(bamfilestring, regexStr=r"\w{2}:\w:"))


@cython.returns(dict)
def GetBamTagTypeDict(cystr bamfile):
    cdef list i
    return {i[0]: i[1] for
            i in [f.split(":") for
                  f in GetBamTagTypes(bamfile)]}

TagTypeDict = {"PV": "B", "AF": "f", "BS": "Z", "FA": "B",
               "FM": "i", "FP": "i", "MQ": "i", "ND": "i",
               "NF": "f", "NM": "i", "RP": "Z", "SC": "Z",
               "SF": "f", "SV": "Z", "X0": "i", "X1": "i",
               "XM": "i", "YA": "Z", "YM": "i", "YO": "Z",
               "YQ": "i", "YR": "i", "YX": "i", "MP": "A",
               "PM": "Z", "MA": "Z", "ot": "i", "mp": "i",
               "om": "i"}


cdef class pFastqFile(object):
    """
    Contains a handle for a kseq.h wrapper and converts each FastqProxy
    to a pFastqProxy
    """

    def __init__(self, object arg):
        if(isinstance(arg, str)):
            self.handle = pysam.cfaidx.FastqFile(arg, persist=False)
        elif(isinstance(arg, pysam.cfaidx.FastqFile)):
            self.handle = pysam.cfaidx.FastqFile(arg.filename, persist=False)
        else:
            raise IllegalArgumentError(
                "pFastqFile can only be initialized from a string "
                "(path to fastq) or a FastqFile object.")

    def __iter__(self):
        return self

    @cython.returns(pFastqProxy_t)
    def __next__(self):
        """
        python version of next().
        """
        return pFastqProxy.fromFastqProxy(self.handle.next())

    def refresh(self):
        self.handle = pysam.FastqFile(self.handle.filename)

    cpdef close(self):
        self.handle.close()


cpdef bint ReadsOverlap(
        AlignedSegment_t read1,
        AlignedSegment_t read2):
    """cpdef wrapper of cReadsOverlap
    """
    return cReadsOverlap(read1, read2)


cdef extern from "math.h":
    double sqrt(double m)


@cython.boundscheck(False)
def cyStdInt(ndarray[np.int64_t, ndim=1] a):
    """Taken from Notes On Cython, with some minor changes.
    """
    return cyOptStdDev_(a.astype(np.float64))


@cython.boundscheck(False)
def cyStdFlt(ndarray[np.float64_t, ndim=1] a):
    """Taken from Notes On Cython, with some minor changes.
    """
    return cyOptStdDev_(a)


@cython.boundscheck(False)
cdef double cyOptStdDev_(ndarray[np.float64_t, ndim=1] a):
    """Taken from Notes On Cython, with some minor changes.
    """
    cdef Py_ssize_t i
    cdef Py_ssize_t n = a.shape[0]
    cdef double m = 0.0
    for i in range(n):
        m += a[i]
    m /= n
    cdef double v = 0.0
    for i in range(n):
        v += (a[i] - m)**2
    return sqrt(v / n)

PhageRefIDDict = {0: 'gi|215104|gb|J02459.1|LAMCG'}
