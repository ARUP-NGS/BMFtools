# cython: c_string_type=str, c_string_encoding=ascii
import abc
from array import array
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
from string import maketrans
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
import multiprocessing as mp
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

@cython.boundscheck(False)
@cython.initializedcheck(False)
@cython.wraparound(False)
@cython.returns(cystr)
def prevcmp(cystr src, uint64_t l):
    cdef char buf[200]
    revcmp(buf, <char *>src, l)
    return <bytes>buf[:l]


@cython.returns(bint)
def lisnone(x):
    return x[1] is None


@cython.returns(tuple)
def lreverse(x):
    return (x[1], x[0])


@cython.returns(int)
def linsertsize(x):
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
            strLen = len(tmpList[tmpInt])
            if strLen < kmerLen:
                tmpList[tmpInt] = "A" * (kmerLen - strLen) + tmpList[tmpInt]
        return tmpList


@memoize
@cython.returns(cystr)
def MemoRevCmp(cystr seq):
    return RevCmp(seq)

DNA_CODON_TABLE = maketrans("ACGTN", "TGCAN")

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

    @classmethod
    def empty(cls):
        return cls("", "", "", "")

    cdef cystr tostring(self):
        return "@%s %s\n%s\n+\n%s\n" % (self.name, self.comment,
                                        self.sequence, self.quality)

    def __str__(self):
        return self.tostring()

    cpdef dict get_tag_dict(self):
        return {el.split(":")[0]: el.split(":")[2] for
                el in self.comment.split("\t")}

    cpdef int get_int_tag(self, cystr key):
        return int(self.get_tag_dict()[key])

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


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef cystr getBS(pFastqProxy_t read):
    return cGetBS(read)


@cython.boundscheck(False)
@cython.wraparound(False)
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


def BwaMemCall(R1, R2, ref="default",
               outBAM="default", path="default",
               bint coorsort=True, bint u=False,
               sortMem="6G", cystr opts=None,
               bint dry_run=False, bint sam=False,
               int threads=4):
    """
    :param: R1 - [cystr/arg] - path to input fastq for read 1
    :param: R2 - [cystr/arg] - path to input fastq for read 2
    :param: ref [cystr/kwarg/"default"] path to reference index base
    :param: outBAM - [cystr/kwarg/"default"] - path to output bam.
    Set to 'stdout' to emit to stdout.
    :param: path - [cystr/kwarg/"default"] - absolute path to bwa executable.
    :param: coorsort [bint/kwarg/True] - whether or not to coordinate sort
    :param: u [bint/kwarg/False] - emit uncompressed bam.
    Override default bwa path (bwa) if necessary.
    :param: sortMem - [cystr/kwarg/"6G"] - sort memory limit for samtools
    :param: opts - [cystr/kwarg/"-t 4 -v 1 -Y -T 0"] - optional arguments
    to provide to bwa for alignment.
    :param: dry_run - [bint/kwarg/False] - flag to return the command string
    rather than calling it.
    :param: threads - [int/kwarg/4]
    :returns: [cystr] - path to outBAM if writing to file, "stdout" if
    emitting to stdout.
    """
    if(opts is None):
        opts = "-t %i -v 1 -Y " % threads
    if(path == "default"):
        path = "bwa"
    if(outBAM == "default"):
        outBAM = ".".join(R1.split(".")[0:-1]) + ".mem.bam"
    if(ref == "default"):
        raise Tim("Reference file index required for alignment!")
    uuidvar = str(uuid.uuid4().get_hex().upper()[0:8])
    opt_concat = ' '.join(opts.split())
    cStr = "%s mem -C %s %s %s %s " % (path, opt_concat, ref, R1, R2)
    if(sam is False):
        if(coorsort):
            compStr = " -l 0 " if(u) else ""
            cStr += " | samtools sort -m %s -O bam -T %s %s -" % (sortMem,
                                                                  uuidvar,
                                                                  compStr)
            if(outBAM not in ["stdout", "-"]):
                cStr += " -o %s" % outBAM
        else:
            cStr += (" | samtools view -Sbhu - " if(
                u) else " | samtools view -Sbh -")
            if(outBAM not in ["stdout", "-"]):
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


def CoorSortAndIndexBam(inBAM, prefix="MetasyntacticVar",
                        outBAM="default",
                        uuid="true", sortMem="4G",
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
                  " -@ {} {} -m {}".format(threads, inBAM, sortMem))
    printlog("About to call sort command: {}".format(CommandStr))
    subprocess.check_call(shlex.split(CommandStr))
    printlog("Now indexing.")
    subprocess.check_call(shlex.split("samtools index {}".format(outBAM)))
    if(delete):
        subprocess.check_call(["rm", inBAM])
    printlog("Output BAM is: " + outBAM)
    return outBAM


def NameSort(inBAM, outBAM="default", prefix="MetasyntacticVar",
             uuid="true", threads="4", sortMem="4G"):
    # If uuid is either a boolean true or is a string containing true,
    # then a random string is generated for the output
    if(str(uuid).lower() == "true"):
        import uuid
        prefix += str(uuid.uuid4().get_hex().upper()[0:8])
    if(outBAM == "default"):
        outBAM = '.'.join(inBAM.split('.')[0:-1]) + '.NameSort.bam'
    CommandStr = ("samtools sort -T {} -O bam -o {}".format(prefix, outBAM) +
                  " -@ {} -m {} -n {}".format(threads, sortMem, inBAM))
    printlog("About to call sort command: {}".format(CommandStr))
    subprocess.check_call(shlex.split(CommandStr))
    printlog("Namesort successful, sorted bam available at: {}".format(outBAM))
    return outBAM


def NameSortAndFixMate(inBAM, outBAM="default", prefix="MetasyntacticVar",
                       uuid="true", threads="4", sortMem="4G",
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


cdef class pPileupRead:
    """
    Python container for the PileupRead proxy in pysam
    """

    cpdef opt(self, cystr arg):
        return self.alignment.opt(arg)

    def __init__(self, pysam.calignedsegment.PileupRead PileupRead):
        cdef py_array BQs, FAs
        self.alignment = PileupRead.alignment
        self.indel = PileupRead.indel
        self.level = PileupRead.level
        self.query_position = PileupRead.query_position
        self.name = self.alignment.qname
        self.BaseCall = self.alignment.seq[self.query_position]
        try:
            BQs = self.alignment.opt("PV")
            self.BQ = BQs[self.alignment.inferred_length -
                          self.query_position - 1]
        except TypeError:
            # Must be old code.
            BQs = array('l', np.array(self.alignment.opt("PV").split(","),
                                      dtype=np.int64))
            self.BQ = BQs[self.alignment.inferred_length -
                          self.query_position - 1]
        except KeyError:
            # Must not have PV tags.
            BQs = self.alignment.query_qualities
            self.BQ = BQs[self.query_position]
        try:
            FAs = self.alignment.opt("FA")
        except TypeError:
            # Must be old code.
            FAs = array('l', np.array(self.alignment.opt("FA").split(","),
                                      dtype=np.int64))
        except KeyError:
            # Must not have FA tags.
            warnings.warn("Missing FA tag in pPileupRead initialization.")
            FAs = array('l', [1] * self.alignment.inferred_length)
        if(self.alignment.is_reverse):
            self.FA = FAs[self.alignment.inferred_length -
                          self.query_position - 1]
        else:
            self.FA = FAs[self.query_position]
        self.MBQ = nmax(BQs)
        self.FM = self.alignment.opt("FM")
        self.MQ = self.alignment.mapping_quality


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


@cython.returns(dict)
def parseConfig(cystr string):
    """
    Parses in a file into a dictionary of key value pairs.
    Key is line.strip().split("=")[0].
    Value is line.strip().split("=")[1].
    Any further values are ignored.
    New with BMFTools v0.0.5.2 (or so?): # comment the rest of a line out.
    """
    return {line.split("=")[0].strip(): line.split("=")[1].strip() for
            line in [l.strip().split("#")[0] for l in
                     open(string, "r").readlines()
                     if l[0] != "#"]}


@cython.returns(dict)
def ReadListToCovCounter(reads, int minClustDepth=3,
                         int minPileupLen=10):
    """
    Makes a Counter object of positions covered by a set of reads.
    Only safe at this point for intrachromosomal rearrangements!
    """
    return cyfreq(reduce(lambda x, y: x + y,
                         [r.get_reference_positions() for r in reads]))


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


def BuildEEModels(cystr f1, cystr f2,
                  np.longdouble_t outliers_fraction=0.1,
                  np.longdouble_t contamination=0.005,
                  int window=20):
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
        cdef cystr key, value
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


def MergeBamList(bamlist, picardpath="default", sortMem="-Xmx6G",
                 outbam="default"):
    """
    Merges a list of BAMs. Used for merging discordant read bams for
    parallelized VCF calls.
    """
    if(isinstance(bamlist, str)):
        bamlist = bamlist.split(":")
    if(outbam == "default"):
        outbam = bamlist[0].split(".")[0:-1] + ".merged.bam"
    commandStr = "java %s -jar %s AS=true" % (sortMem, picardpath)
    commandStr += " I=" + " I=".join(bamlist)
    commandStr += " O=%s" % outbam
    check_call(shlex.split(commandStr))
    return outbam


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


cdef inline bint TestBarcode(char *BS, int8_t hpLimit, int bLen) nogil:
    cdef char nuc, last
    cdef int run, index
    last = 0
    run = 0
    for index in range(bLen):
        nuc = BS[index]
        if nuc == 78:
            return False
        if nuc == last:
            run += 1
        else:
            run = 0
    return run < hpLimit
