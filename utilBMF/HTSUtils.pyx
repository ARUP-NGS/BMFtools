# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True
from __future__ import division
import abc
from Bio.Seq import Seq
from copy import copy as ccopy
from cytoolz import map as cmap, memoize, frequencies as cyfreq
from functools import partial
from itertools import groupby, tee, chain, combinations, product
from MawCluster.Probability import GetCeiling
from numpy import any as npany
from numpy import concatenate as nconcatenate
from numpy import less as nless
from numpy import max as nmax
from numpy import mean as nmean
from numpy import min as nmin
from numpy import std as nstd
from operator import iadd as oia
from operator import itemgetter as oig
from pysam.calignmentfile import AlignedSegment as pAlignedSegment
from re import compile as regexcompile
from subprocess import check_output, check_call, CalledProcessError
from .ErrorHandling import *
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

cimport numpy as np
cimport pysam.cfaidx
cimport pysam.calignmentfile
cimport pysam.TabProxies
ctypedef pysam.calignmentfile.AlignedSegment cAlignedSegment
ctypedef pysam.calignmentfile.PileupRead cPileupRead
ctypedef Insertion Insertion_t
ctypedef Deletion Deletion_t
ctypedef AbstractIndelContainer AbstractIndelContainer_t
ctypedef IndelQuiver IndelQuiver_t
oig1 = oig(1)
oig0 = oig(0)
cfi = chain.from_iterable


def l1(x):
    return x[1]


def lisnone(x):
    return x[1] is None


def lreverse(x):
    return (x[1], x[0])


def linsertsize(x):
    return x.insert_size


def printlog(string, level=logging.INFO):
    Logger = logging.getLogger("Primarylogger")
    if(level == logging.DEBUG):
        Logger.debug(string.replace(
            "\t", "\\t").replace("\n", "\\n").replace(
                "'", "\'").replace('"', '\\"'))
        return
        # Doesn't print to string if set to debug mode.
    elif(level == logging.INFO):
        Logger.info(string.replace(
            "\t", "\\t").replace("\n", "\\n").replace(
                "'", "\'").replace('"', '\\"'))
    elif(level == logging.WARNING):
        Logger.warning(string.replace(
            "\t", "\\t").replace("\n", "\\n").replace(
                "'", "\'").replace('"', '\\"'))
    else:
        Logger.critical(string.replace(
            "\t", "\\t").replace("\n", "\\n").replace(
            "'", "\'").replace('"', '\\"'))
    sys.stderr.write(string.replace(
        "\t", "\\t").replace("\n", "\\n").replace(
        "'", "\'").replace('"', '\\"') + "\n")
    return

pl = printlog

CmpDict = {"A": "T", "C": "G", "G": "C", "T": "A"}


@cython.returns(cython.str)
def RevCmp(cython.str seq):
    return "".join([CmpDict[i] for i in list(seq)])[::-1]


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

ChrToPysamDict = {}
for i in xrange(22):
    ChrToPysamDict[str(i + 1)] = i
ChrToPysamDict["*"] = -1
ChrToPysamDict["X"] = 22
ChrToPysamDict["Y"] = 23
ChrToPysamDict["MT"] = 24
ChrToPysamDict["GL000207.1"] = 25
ChrToPysamDict["GL000226.1"] = 26
ChrToPysamDict["GL000229.1"] = 27
ChrToPysamDict["GL000231.1"] = 28
ChrToPysamDict["GL000210.1"] = 29
ChrToPysamDict["GL000239.2"] = 30
ChrToPysamDict["GL000235.1"] = 31
ChrToPysamDict["GL000201.1"] = 32
ChrToPysamDict["GL000247.1"] = 33
ChrToPysamDict["GL000245.1"] = 34
ChrToPysamDict["GL000197.1"] = 35
ChrToPysamDict["GL000203.1"] = 36
ChrToPysamDict["GL000246.1"] = 37
ChrToPysamDict["GL000249.1"] = 38
ChrToPysamDict["GL000196.1"] = 39
ChrToPysamDict["GL000248.1"] = 40
ChrToPysamDict["GL000244.1"] = 41
ChrToPysamDict["GL000238.1"] = 42
ChrToPysamDict["GL000202.1"] = 43
ChrToPysamDict["GL000234.1"] = 44
ChrToPysamDict["GL000232.1"] = 45
ChrToPysamDict["GL000206.1"] = 46
ChrToPysamDict["GL000240.1"] = 47
ChrToPysamDict["GL000236.1"] = 48
ChrToPysamDict["GL000241.1"] = 49
ChrToPysamDict["GL000243.1"] = 50
ChrToPysamDict["GL000242.1"] = 51
ChrToPysamDict["GL000230.1"] = 52
ChrToPysamDict["GL000237.1"] = 53
ChrToPysamDict["GL000233.1"] = 54
ChrToPysamDict["GL000204.1"] = 55
ChrToPysamDict["GL000198.1"] = 56
ChrToPysamDict["GL000208.1"] = 57
ChrToPysamDict["GL000191.1"] = 58
ChrToPysamDict["GL000227.1"] = 59
ChrToPysamDict["GL000228.1"] = 60
ChrToPysamDict["GL000214.1"] = 61
ChrToPysamDict["GL000221.1"] = 62
ChrToPysamDict["GL000209.1"] = 63
ChrToPysamDict["GL000218.1"] = 64
ChrToPysamDict["GL000220.1"] = 65
ChrToPysamDict["GL000213.1"] = 66
ChrToPysamDict["GL000211.1"] = 67
ChrToPysamDict["GL000199.1"] = 68
ChrToPysamDict["GL000217.1"] = 69
ChrToPysamDict["GL000216.1"] = 70
ChrToPysamDict["GL000215.1"] = 71
ChrToPysamDict["GL000205.1"] = 72
ChrToPysamDict["GL000219.1"] = 73
ChrToPysamDict["GL000224.1"] = 74
ChrToPysamDict["GL000223.1"] = 75
ChrToPysamDict["GL000195.1"] = 76
ChrToPysamDict["GL000212.1"] = 77
ChrToPysamDict["GL000222.1"] = 78
ChrToPysamDict["GL000200.1"] = 79
ChrToPysamDict["GL000193.1"] = 80
ChrToPysamDict["GL000194.1"] = 81
ChrToPysamDict["GL000225.1"] = 82
ChrToPysamDict["GL000192.1"] = 83


@cython.returns(dict)
def GetPysamToChrDict(cython.str alignmentFileText):
    """
    Input variable contains the "text" attribute from a pysam.AlignmentFile
    object.
    """
    global PysamToChrDict
    PysamToChrDict = dict(list(enumerate(
        [i.replace("SN:", "").split("\t")[1] for i in
         alignmentFileText.split('\n') if i[0:3] == "@SQ"])))
    global ChrToPysamDict
    ChrToPysamDict = {PysamToChrDict[key]: key for key in
                      PysamToChrDict.iterkeys()}
    return


@cython.returns(dict)
def GetPysamToChrDictFromAlignmentFile(
        pysam.calignmentfile.AlignmentFile alignmentfileObj):
    """
    Returns a dictionary of pysam reference numbers to contig names.
    """
    return dict(list(enumerate(alignmentfileObj.references)))


@cython.returns(dict)
def GetChrToPysamDictFromAlignmentFile(alignmentfileObj):
    """
    Returns a dictionary of contig names to pysam reference numbers.
    """
    assert isinstance(alignmentfileObj, pysam.calignmentfile.AlignmentFile)
    return dict(list(cmap(lreverse,
                     list(enumerate(alignmentfileObj.references)))))


@cython.returns(dict)
def GetBidirectionalPysamChrDict(alignmentfileObj):
    """
    Returns a dictionary of contig names to pysam reference numbers
    and vice versa - bi-directional.
    """
    assert isinstance(alignmentfileObj, pysam.calignmentfile.AlignmentFile)
    refList = list(enumerate(alignmentfileObj.references))
    return dict(list(cmap(lreverse, refList) + refList))


class pFastqProxy:
    """
    Python container for pysam.cfaidx.FastqProxy with persistence.
    """
    @cython.locals(FastqProxyObj=pysam.cfaidx.FastqProxy)
    def __init__(self, FastqProxyObj):
        self.comment = FastqProxyObj.comment
        self.quality = FastqProxyObj.quality
        self.sequence = FastqProxyObj.sequence
        self.name = FastqProxyObj.name

    def __str__(self):
        return ("@" + self.name + " " + self.comment +
                "\n" + self.sequence + "\n+\n" + self.quality + "\n")


def FacePalm(string):
    Str = ("............................................________ "
           "\n....................................,.-'\"...................`"
           "`~. \n.............................,.-\"........................"
           "...........\"-., \n.........................,/.................."
           ".............................\":, \n.....................,?....."
           "................................................., \n..........."
           "......../........................................................"
           "...,} \n................./......................................."
           "...............,:`^`..} \n.............../......................."
           "............................,:\"........./ \n..............?....."
           "__.........................................:`.........../ \n....."
           "......../__.(.....\"~-,_..............................,:`........"
           "../ \n.........../(_....\"~,_........\"~,_....................,:`"
           "........_/ \n..........{.._$;_......\"=,_.......\"-,_.......,.-~-"
           ",},.~\";/....} \n...........((.....*~_.......\"=-._......\";,,./`"
           "..../\"............../ \n...,,,___.`~,......\"~.,................"
           "....`.....}............../ \n............(....`=-,,.......`......"
           "..................(......;_,,-\" \n............/.`~,......`-....."
           "................................/ \n.............`~.*-,.........."
           "...........................|,./.....,__ \n,,_..........}.>-._...."
           "...............................|..............`=~-, \n.....`=~-,_"
           "_......`,................................. \n...................`"
           "=~-,,.,............................... \n........................"
           "........`:,,...........................`..............__ \n......"
           "...............................`=-,...................,%`>--==`` "
           "\n........................................_..........._,-%.......`"
           " \n...................................,")
    print(Str)
    if(isinstance(string, str)):
        raise ThisIsMadness(string)
    raise ThisIsMadness("WHAT YOU SAY")


@cython.returns(cython.bint)
def is_read_softclipped(read):
    """
    Simply returns whether or not a read is soft-clipped
    """
    if(read.cigarstring is None):
        return False
    if("S" in read.cigarstring):
        return True
    return False


@cython.returns(cython.bint)
@cython.locals(minLen=cython.long)
def ReadPairIsDuplex(readPair, minShare="default"):
    """
    If minShare is an integer, require that many nucleotides
    overlapping to count it as duplex.
    Defaults to sharing at least half.
    """
    if(readPair.read1_contig != readPair.read2_contig):
        return False
    if(isinstance(minShare, int)):
        minLen = minShare
    elif(isinstance(minShare, float)):
        minLen = int(minShare * readPair.read1.query_length)
    elif(minShare == "default"):
        minLen = readPair.read1.query_length // 2
    else:
        raise ThisIsMadness("minShare parameter required. Integer for "
                            "an absolute number of bases overlapped re"
                            "quired, float for a fraction of read length.")
    return sum([x == 2 for x in
                cyfreq(readPair.read1.get_reference_positions() +
                       readPair.read2.get_reference_positions()).values()]
               ) >= minLen


def BwaswCall(fq1, fq2, ref="default", outBAM="default"):
    if(ref == "default"):
        raise ThisIsMadness("ref required to call bwasw.")
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


def align_bwa_aln_addRG(R1, R2, ref="default", opts="", outBAM="default",
                        picardPath="default", RG="default",
                        PL="ILLUMINA", SM="default", ID="default",
                        CN="default"):
    """
    Aligns a set of paired-end reads using bwa aln. Defaults to 4 threads.
    In order to make BAMs compatible with both GATK and pysam,

    """
    if(ref == "default"):
        FacePalm("Reference file index required for alignment!")
    if(picardPath == "default"):
        FacePalm("Picard jar path required for adding read groups!")
    if(opts == ""):
        opts = "-n 3 -t 4"
    if(outBAM == "default"):
        outBAM = '.'.join(R1.split('.')[0:-1]) + ".aln.bam"
    outSAM = outBAM.replace("bam", "sam")
    # Note: ID has to be "bwa" so that it passes the SAM validation that
    # the hsdjdk has, which is required for either GATK realignment or
    # ABRA realignment.
    str(uuid.uuid4().get_hex().upper()[0:8])
    R1Sai = R1 + ".tmp.sai"
    R2Sai = R2 + ".tmp.sai"
    alnStr1 = ("bwa aln " + opts + " " + " ".join([ref, R1]) +
               " > " + R1Sai)
    PipedShellCall(alnStr1)
    alnStr2 = ("bwa aln " + opts + " " + " ".join([ref, R2]) +
               " > " + R2Sai)
    PipedShellCall(alnStr2)
    sampeStr = ("bwa sampe " + " ".join([ref, R1Sai, R2Sai, R1, R2]) +
                "  > " + outSAM)
    printlog("bwa aln string: {}".format(sampeStr))
    PipedShellCall(sampeStr)
    AddReadGroupsPicard(outSAM, outBAM=outBAM, SM=SM, ID=ID, PL=PL,
                        CN=CN, picardPath=picardPath)
    os.remove(R1Sai)
    os.remove(R2Sai)
    return outBAM


def align_bwa_aln(R1, R2, ref="default", opts="", outBAM="default"):
    """
    Aligns a set of paired-end reads using bwa aln. Defaults to 4 threads.
    In order to make BAMs compatible with both GATK and pysam,

    """
    if(ref == "default"):
        FacePalm("Reference file index required for alignment!")
    if(opts == ""):
        opts = "-n 3 -t 4"
    if(outBAM == "default"):
        outBAM = '.'.join(R1.split('.')[0:-1]) + ".aln.bam"
    outSAM = outBAM.replace("bam", "sam")
    # Note: ID has to be "bwa" so that it passes the SAM validation that
    # the hsdjdk has, which is required for either GATK realignment or
    # ABRA realignment.
    str(uuid.uuid4().get_hex().upper()[0:8])
    R1Sai = R1 + ".tmp.sai"
    R2Sai = R2 + ".tmp.sai"
    alnStr1 = ("bwa aln " + opts + " " + " ".join([ref, R1]) +
               " > " + R1Sai)
    PipedShellCall(alnStr1)
    alnStr2 = ("bwa aln " + opts + " " + " ".join([ref, R2]) +
               " > " + R2Sai)
    PipedShellCall(alnStr2)
    sampeStr = ("bwa sampe " + " ".join([ref, R1Sai, R2Sai, R1, R2]) +
                " | samtools view -h - > " + outBAM)
    printlog("bwa aln string: {}".format(sampeStr))
    PipedShellCall(sampeStr)
    os.remove(R1Sai)
    os.remove(R2Sai)
    return outBAM


def align_bwa_mem_addRG(R1, R2, ref="default", opts="", outBAM="default",
                        path="default", picardPath="default",
                        PL="ILLUMINA", SM="default", ID="default",
                        CN="default"):
    """
    Aligns a set of paired-end
    reads to a reference
    with provided options using bwa mem.
    Defaults to 4 threads, silent alignment, listing
    supplementary alignments, and
    writing each reads' alignment,
    regardless of mapping quality.
    """
    if(opts == ""):
        opts = '-t 4 -v 1 -Y -M -T 0'
    if(outBAM == "default"):
        outBAM = ".".join(R1.split(".")[0:-1]) + ".mem.bam"
    outSAM = outBAM.replace(".bam", ".sam")
    if(ref == "default"):
        FacePalm("Reference file index required for alignment!")
    if(picardPath == "default"):
        FacePalm("Path to picard jar required for adding RG!")
    opt_concat = ' '.join(opts.split())
    RGString = "@RG\tID:bwa SM:%s" % SM
    if(path == "default"):
        command_str = ("bwa mem %s %s %s %s " % (opt_concat, ref, R1, R2) +
                       " > {}".format(outSAM))
    else:
        command_str = (path + " mem %s %s %s " % (opt_concat, ref, R1) +
                       "{} > {}".format(R2, outSAM))
    # command_list = command_str.split(' ')
    printlog(command_str)
    PipedShellCall(command_str)
    AddReadGroupsPicard(outSAM, outBAM=outBAM, SM=SM, ID=ID, PL=PL,
                        CN=CN, picardPath=picardPath)
    return outBAM


def align_bwa_mem(R1, R2, ref="default", opts="", outBAM="default",
                  path="default"):
    """
    Aligns a set of paired-end
    reads to a reference
    with provided options using bwa mem.
    Defaults to 4 threads, silent alignment, listing
    supplementary alignments, and
    writing each reads' alignment,
    regardless of mapping quality.
    """
    if(opts == ""):
        opts = '-t 4 -v 1 -Y -M -T 0'
    if(outBAM == "default"):
        outBAM = ".".join(R1.split(".")[0:-1]) + ".mem.bam"
    if(ref == "default"):
        FacePalm("Reference file index required for alignment!")
    opt_concat = ' '.join(opts.split())
    if(path == "default"):
        command_str = ("bwa mem %s %s %s %s " % (opt_concat, ref, R1, R2) +
                       " | samtools view -Sbh - > {}".format(outBAM))
    else:
        command_str = (path + " mem %s %s %s " % (opt_concat, ref, R1) +
                       "{} | samtools view -Sbh - > {}".format(R2, outBAM))
    # command_list = command_str.split(' ')
    printlog(command_str)
    PipedShellCall(command_str, delete=True)
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
    PipedShellCall(command_str, delete=True)
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


def PipedShellCall(commandStr, delete=True, silent=False):
    PipedShellCallFilename = "PipedShellCall{}.sh".format(
        str(uuid.uuid4().get_hex().upper()[0:8]))
    if silent is False:
        printlog("Command string: {}".format(commandStr), level=logging.DEBUG)
    open(PipedShellCallFilename, "w").write(commandStr)
    subprocess.check_call(['bash', PipedShellCallFilename])
    if(delete):
        try:
            os.remove(PipedShellCallFilename)
        except OSError:
            printlog("Piped shell call deletion failed. Oh well!",
                     level=logging.DEBUG)
    return commandStr


def CustomRefBowtiePaired(mergedFq,
                          ref,
                          output="default",
                          barLen="default",
                          bowtiePath="bowtie",
                          mismatchLimit="default"):
    if(output == "default"):
        output = mergedFq.split('.')[0] + '.mergingFamilies.sam'
    if(barLen == "default"):
        FacePalm("Barcode length must be set. Abort mission!")
    if(bowtiePath == "bowtie"):
        printlog("Defaulting to bowtie for path to executable.")
    elif("2" in bowtiePath.split("/")[-1]):
        FacePalm("Do not use bowtie2!")
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
            raise ThisIsMadness("OMGZ")
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
        FacePalm("VCFLineContainedInBed requires a VCFRecord object!")
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
        raise ThisIsMadness("A reference file path must be provided!")
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
        raise ThisIsMadness("A reference file path must be provided!")
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


@cython.locals(sortAndIndex=cython.bint)
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
        outBAM = ".".join(bamlist[0].split(".")[-1]) + ".merged.bam"
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

    def __init__(self, pysam.calignmentfile.AlignedSegment read1,
                 pysam.calignmentfile.AlignedSegment read2):
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
        # TODO: write a script to create an array of soft-clipped sequences
        # from each read

    @cython.returns(cython.long)
    def NumOverlappingBed(self, list bedLines=[]):
        try:
            assert isinstance(bedLines[0], str) and isinstance(
                bedLines[1], int)
        except AssertionError:
            FacePalm("Sorry, bedLines must be in ParseBed format.")
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


cdef class pPileupRead:
    """
    Python container for the PileupRead proxy in pysam
    """

    def __init__(self, pysam.calignmentfile.PileupRead PileupRead):
        self.alignment = PileupRead.alignment
        self.indel = PileupRead.indel
        self.level = PileupRead.level
        self.query_position = PileupRead.query_position
        self.name = self.alignment.qname
        self.BaseCall = self.alignment.seq[self.query_position]


cdef class PileupReadPair:

    """
    Holds both bam record objects in a pair of pileup reads.
    Currently, one read unmapped and one read soft-clipped are
    both marked as soft-clipped reads.
    Accepts a list of length two as input.
    """

    def __cinit__(self, tuple readlist):
        cdef pPileupRead_t read1
        cdef pPileupRead_t read2
        read1, read2 = readlist[0], readlist[1]
        try:
            assert len(readlist) == 2
        except AssertionError:
            pl("repr(readlist): %s" % repr(readlist))
            raise ThisIsMadness(
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


def GetReadPair(inHandle):
    """
    Simply contains both pairs of reads in an object
    """
    read1 = inHandle.next()
    read2 = inHandle.next()
    try:
        assert read1.query_name == read2.query_name
    except AssertionError:
        FacePalm("These two reads have "
                 "different query names. Abort!")
    return ReadPair(read1, read2)


def ReadsOverlap(read1, read2):
    if(read1.reference_id != read2.reference_id):
        return False
    if(read1.reference_start > read2.reference_end or
       read1.reference_end < read2.reference_start):
        return False
    else:
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
            FacePalm("Minimum family size must be castable to int!")
        RecordsArray = [rec for rec in RecordsArray if
                        rec.opt("FM") >= minFamSize]
    inHandle.close()
    return RecordsArray


def LoadReadPairsFromFile(inBAM, SVTag="default",
                          minMQ=0, minBQ=0):
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
        return sorted(RecordsArray, key=linsertize)
    else:
        return RecordsArray


def WritePairToHandle(ReadPair, handle="default"):
    """
    Writes a pair to a file handle.
    """
    assert isinstance(handle, pysam.calignmentfile.AlignmentFile)
    handle.write(ReadPair.read1)
    handle.write(ReadPair.read2)
    return True


@cython.returns(list)
def ParseBed(cython.str bedfile):
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
def parseConfig(cython.str string):
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
def ReadListToCovCounter(reads, cython.long minClustDepth=3,
                         cython.long minPileupLen=10):
    """
    Makes a Counter object of positions covered by a set of reads.
    Only safe at this point for intrachromosomal rearrangements!
    """
    return cyfreq(reduce(lambda x, y: x + y,
                         [r.get_reference_positions() for r in reads]))


@cython.returns(dict)
def ReadPairListToCovCounter(list ReadPairList, cython.long minClustDepth=5,
                             cython.long minPileupLen=10):
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
            raise ThisIsMadness("Soft-clipped seq is ambiguous "
                                "without a contig.")
        if isinstance(is_reverse, bool) is False:
            raise ThisIsMadness("SoftClippedSeq needs to know to which"
                                "strand it is mapped.")
        self.seq = seq
        self.contig = contig
        self.is_reverse = is_reverse

    def RCChangeStrand(self):
        self.seq = str(Seq(self.seq).reverse_complement())
        self.is_reverse = not self.is_reverse


def SplitSCRead(pysam.calignmentfile.AlignedSegment read):
    try:
        assert "S" in read.cigarstring
    except AssertionError:
        FacePalm("You can't split a read by soft-clipping "
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


@cython.returns(list)
def CreateIntervalsFromCounter(dict CounterObj, cython.long minPileupLen=0,
                               cython.str contig="default",
                               bedIntervals="default",
                               cython.long mergeDist=0,
                               cython.long minClustDepth=5):
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
        FacePalm("contig required for this function!")
    for k, g in groupby(
            enumerate(sorted(CounterObj.iterkeys())), lix):
        posList = list(cmap(oig1, g))
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

# Storing the numconv functions for easy application.
Base64ToInt = numconv.NumConv(64).str2int
Int2Base64 = numconv.NumConv(64).int2str

ph2chrDict = {}
for i in xrange(100):
    ph2chrDict[i] = chr(i + 33)
chr2ph = {ph2chrDict[key]: key for key in
          ph2chrDict.iterkeys()}

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


@cython.locals(rq1=cython.int, n=cython.int)
def CigarToQueryIndices(cigar):
    """
    Returns a list of lists of positions within the read.
    sublist[0] is the cigarOp, sublist[1] is the list of positions.
    """
    if(cigar is None):
        raise ThisIsMadness("I can't get query indices "
                            "if there's no cigar string.")
    try:
        assert isinstance(cigar[0], tuple)
    except AssertionError:
        pl("Invalid argument - cigars are lists of tuples, but the first item"
           " in this list is not a tuple!")
    tuples = []
    c = list(cmap(oig1, cigar))
    cumSum = [sum(c[:i + 1]) for i in xrange(len(c))]
    for n, entry in enumerate(cigar):
        if n == 0:
            tuples.append((entry[0], range(entry[1])))
        else:
            tuples.append((entry[0], range(cumSum[n - 1], cumSum[n])))
    return tuples


def GetQueryIndexForCigarOperation(pysam.calignmentfile.AlignedSegment read,
                                   cython.long cigarOp=-1):
    """
    Returns a list of lists of positions within each read
    which match the given cigarOp. The cigarOp must be an integer.
    See pysam's documentation for details.
    """
    if(read.cigar is None):
        raise ThisIsMadness("I can't get query indices "
                            "if there's no cigar string.")
    try:
        assert(cigarOp in range(9))
    except AssertionError:
        raise ThisIsMadness("Please read the doc string for "
                            "GetQueryIndexForCigarOperation. Invalid cigarOp")
    assert isinstance(read, pysam.calignmentfile.AlignedSegment)
    # Get positions in read which match cigar operation.
    QueryCigar = CigarToQueryIndices(read.cigar)
    filtCigar = [i for i in QueryCigar if i[0] == cigarOp]
    return filtCigar


def GetReadNucsForFiltCigar(filtCigar,
                            pysam.calignmentfile.AlignedSegment read):
    """
    Returns a list of strings for nucleotides matching a filtCigar.
    """
    seq = read.seq
    return ["".join([seq[i] for i in g[1]]) for g in filtCigar]


def GetGenomicCoordsForFiltCigar(filtCigar,
                                 pysam.calignmentfile.AlignedSegment read):
    """
    Returns a list of lists for genomic coordinates matching a filtCigar.
    """
    cdef dgap
    dgap = dict(read.get_aligned_pairs())
    return [[dgap[i] for i in g[1]] for g in filtCigar]


def GetGenomicCoordToNucleotideMapForFiltCigar(
        pysam.calignmentfile.AlignedSegment read,
        filtCigar="default"):
    """
    Returns a dictionary of genomic positions: nucleotides for a "filtered"
    cigar and a read. This is used to compare whether or not two reads with
    an insertion agree on the inserted nucleotides.
    """
    assert not isinstance(filtCigar, str)
    seq = read.seq
    dgap = dict(read.get_aligned_pairs())
    try:
        dgapList = [[dgap[i] for i in g[1]] for g in filtCigar]
        seqList = [[seq[i] for i in g[1]] for g in filtCigar]
        return {k: l for k, l in zip(
            [[dgap[i] for i in g[1]] for g in filtCigar],
            [[seq[i] for i in g[1]] for g in filtCigar])}
    except TypeError:
        print(repr(dgap))
        print(repr(filtCigar))
        print(repr(seq))
        if("dgapList" in locals()):
            print(repr(dgapList))
        if("seqList" in locals()):
            print(repr(seqList))
        print(read.cigarstring)
        raise TypeError


@cython.locals(cigarOp=cython.long)
def GetGC2NMapForRead(pysam.calignmentfile.AlignedSegment read,
                      cigarOp=-1):
    """
    Returns a dictionary of genomic positions: nucleotides for a provided
    cigar operation.
    """
    filtCigar = GetQueryIndexForCigarOperation(read, cigarOp=cigarOp)
    return GetGenomicCoordToNucleotideMapForFiltCigar(read,
                                                      filtCigar=filtCigar)


@cython.locals(cigarOp=cython.long)
def GetReadSequenceForCigarOp(pysam.calignmentfile.AlignedSegment read,
                              cigarOp=-1):
    if(read.cigar is None):
        raise ThisIsMadness("Cigar must not be None to get "
                            "sequence for an operation!")
    try:
        assert(cigarOp in range(9))
    except AssertionError:
        raise ThisIsMadness("Please read the doc string for "
                            "GetReadSequenceForCigarOp. Invalid cigarOp")
    filtCigar = GetQueryIndexForCigarOperation(read)
    seqTups = []
    for f in filtCigar:
        if(len(f) == 0):
            continue
        else:
            pass
    raise ThisIsMadness("I haven't finished writing this function. Oops!")


@cython.returns(list)
def GetDeletedCoordinates(pysam.calignmentfile.AlignedSegment read):
    """
    Returns a list of integers of genomic coordinates for deleted bases.
    """
    cdef cython.long k
    assert isinstance(read, pysam.calignmentfile.AlignedSegment)
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
def GetInsertedNucleotides(pysam.calignmentfile.AlignedSegment read):
    """
    """
    cdef cython.long start, end
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
def GetInsertedStrs(pysam.calignmentfile.AlignedSegment read):
    """
    Returns a list of tuples for strings, along with preceding reference base
    and successive reference base position. We can then directly compare these
    tuples between read 1 and read 2 if they cover the same positions.
    Used to determine whether or not DSI is appropriate.
    """
    cdef dict readPosToAlignedPosDict
    cdef cython.long l
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


@cython.returns(dtype128_t)
def FractionSoftClipped(cAlignedSegment read):
    """
    Returns the fraction of a read aligned.
    """
    if(read.cigarstring is None):
        return 0.
    return FractionSoftClippedCigar(read.cigar)


@cython.returns(dtype128_t)
def FractionSoftClippedCigar(list cigar):
    """
    Returns the fraction softclipped directly from tuple
    """
    cdef tuple i
    return 1. * sum(i[1] for i in cigar if i[0] == 4) / sum([i[1] for i
                                                             in cigar])


@cython.returns(dtype128_t)
def FractionAlignedCigar(list cigar):
    """
    Returns the fraction aligned directly from tuple
    """
    cdef tuple i
    return 1. * sum(i[1] for i in cigar if i[0] == 0) / sum([i[1] for i
                                                             in cigar])


@cython.returns(dtype128_t)
def FractionAligned(cAlignedSegment read):
    """
    Returns the fraction of a read aligned.
    """
    if(read.cigarstring is None):
        return 0.
    return FractionAlignedCigar(read.cigar)


def AddReadGroupsPicard(inBAM, RG="default", SM="default",
                        PL="ILLUMINA", CN="default", picardPath="default",
                        outBAM="default", ID="default", LB="default",
                        PU="default"):
    if(picardPath == "default"):
        raise ThisIsMadness("picardPath required to call PicardTools!")
    if(outBAM == "default"):
        outBAM = ".".join(inBAM.split(".")[:-1] + ["addRG", "bam"])
    commandStr = ("java -jar %s AddOrReplaceReadGroups I=" % picardPath +
                  "%s O=%s VALIDATION_STRINGENCY=SILENT " % (inBAM, outBAM) +
                  " CN=%s PL=%s SM=%s ID=%s LB=%s PU=%s" % (CN, PL,
                                                            SM, ID, LB,
                                                            PU))
    printlog("AddReadGroupsPicard commandStr: %s" % commandStr)
    subprocess.check_call(shlex.split(commandStr))
    return outBAM


@cython.locals(outliers_fraction=dtype128_t, contamination=dtype128_t,
               window=cython.long)
def BuildEEModels(f1, f2, outliers_fraction=0.1, contamination=0.005,
                  window=20):
    cdef np.ndarray[dtype128_t, ndim = 1] GAFreqNP = f1
    cdef np.ndarray[dtype128_t, ndim = 1] CTFreqNP = f2
    cdef np.ndarray[dtype128_t, ndim = 1] FreqArray = nconcatenate(GAFreqNP,
                                                                   CTFreqNP)
    ee1 = EllipticEnvelope(contamination=contamination, assume_centered=False)
    ee2 = EllipticEnvelope(contamination=contamination, assume_centered=False)
    GAClassifier = ee1.fit(GAFreqNP)
    CTClassifier = ee2.fit(CTFreqNP)
    pass


def PlotNormalFit(array, outfile="default", maxFreq=0.2):
    import matplotlib as mpl
    mpl.use('Agg')      # With this line = figure disappears
    import matplotlib.pyplot as plt
    import matplotlib.mlab as mlab
    array = array.reshape(-1, 1)
    mu, sigma = nmean(array), nstd(array)
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
    cdef pysam.cfaidx.FastqProxy read
    cdef cython.long famS
    cdef cython.long sumFam
    cdef cython.long numFam
    cdef cython.long numSing
    cdef cython.long sumAll
    cdef dtype128_t meanFamAll
    cdef dtype128_t meanRealFam
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


@cython.locals(n=cython.long)
@cython.returns(list)
def bitfield(n):
    """
    Parses a bitwise flag into an array of 0s and 1s.
    No need - use the & or | tools for working with bitwise flags.
    """
    return [1 if digit == '1' else 0 for digit in bin(n)[2:]]


@cython.returns(cython.str)
def ASToFastqSingle(pysam.calignmentfile.AlignedSegment read):
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
@cython.returns(cython.str)
def ASToFastqPaired(pysam.calignmentfile.AlignedSegment read,
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
def SWRealignAS(pysam.calignmentfile.AlignedSegment read,
                pysam.calignmentfile.AlignmentFile alignmentfileObj,
                cython.str extraOpts="",
                cython.str ref="default", cython.float minAF=0.5):
    """
    Passes the sequence and qualities of the provided read to bwa aln with
    an AlignedSegment as input. Updates the AlignedSegment's fields in-place
    if the result is good.
    Ideally, these records are accessed through a name-sorted AlignmentFile,
    as it throws an index error.
    Might fix this later, but I don't like it very much. If there were cython
    bindings to bwa, that would be the perfect use of it. Ehhh..
    """
    cdef cython.long lAlignedArr, lbf
    cdef cython.float af
    cdef cython.str FastqStr, commandStr
    cdef list alignedArr, letters, numbers, tags, tag
    if(read.opt("AF") == 1.):
        # Nothing wrong a fully-aligned read.
        return read
    try:
        assert ref != "default"
    except AssertionError:
        raise ThisIsMadness("Reference required for AlnRealignAS!")
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
    ra = regexcompile("[A-Z]")
    rn = regexcompile("[0-9]")
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
def makeinfodict(pysam.TabProxies.VCFProxy rec):
    """
    Returns a dictionary of info fields for a tabix VCF Proxy
    """
    return dict([i.split("=") for i in rec.info.split(";")])


@cython.returns(dict)
def makeformatdict(pysam.TabProxies.VCFProxy rec):
    """
    Returns a dictionary of format fields for a tabix VCF Proxy
    """
    return dict(zip(rec.format.split(":"), rec[0].split(":")))


@cython.returns(cython.bint)
def DeaminationConfTest(pysam.TabProxies.VCFProxy rec,
                        dtype128_t ctfreq=-1., dtype128_t conf=1e-3):
    """
    Tests whether or not a given record should pass.
    """
    cdef dict InfoDict
    cdef dict Counts
    cdef cython.long ACR
    cdef cython.long ceiling
    cdef cython.long AC
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
        raise ThisIsMadness("kwarg ctfreq required for DeaminationConfTest!")
    ceiling = GetCeiling(ACR, p=ctfreq, pVal=conf)
    if(AC > ceiling):
        return True
    return False


def PartialDeaminationConfTest(dtype128_t ctfreq, dtype128_t conf=1e-3):
    """
    Returns a function that can be easily used by AbstractVCFProxyFilter.
    """
    return partial(DeaminationConfTest, ctfreq=ctfreq, conf=conf)


class AbstractVCFProxyFilter(object):
    """
    Abstract class which serves as a template for VCF record post-filtering.
    """
    def __init__(self, filterStr, func=FacePalm,
                 key="default", value="*"):
        if(func == FacePalm):
            func("func must be overridden for "
                 "AbstractVCFProxyFilter to work!")
        if(filterStr == "default"):
            FacePalm("filterStr must be set to append or replace the "
                     "filter field.")
        if(hasattr(func, "__call__")):
            self.func = func
        else:
            raise AttributeError("func must be callable!")
        self.filterStr = filterStr
        self.key = key
        self.value = value

    @cython.returns(pysam.TabProxies.VCFProxy)
    def filter(self, pysam.TabProxies.VCFProxy rec):
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


def MakeVCFProxyDeaminationFilter(dtype128_t ctfreq, dtype128_t conf=1e-3,
                                  key="default", value="*"):
    """
    Returns the VCFProxyFilter object I wanted.
    """
    return AbstractVCFProxyFilter("DeaminationNoise",
                                  func=PartialDeaminationConfTest(ctfreq,
                                                                  conf=conf),
                                  key=key, value=value)


def SortBgzipAndTabixVCF(inVCF, outVCF="default",
                         vcflib=True):
    """
    Sorts and tabix indexes a VCF.
    If vcflib is True, use vcfstreamsort from ekg's excellent vcf API.
    """
    if(outVCF == "default"):
        outVCF = ".".join([inVCF.split(".")[-1]]) + ".sort.vcf"
    if(vcflib is False):
        check_call("zcat %s | head -n 1000 | grep '^#' > %s" % (inVCF, outVCF),
                   shell=True)
        check_call("zcat %s | grep -v '^#' >> %s" % (inVCF, outVCF),
                   shell=True)
    else:
        check_call("vcfstreamsort -w 100000 %s > %s" % (inVCF, outVCF))
    check_call(["bgzip", outVCF])
    check_call(["tabix", outVCF + ".gz"])
    return outVCF + ".gz"


@cython.returns(cython.str)
def SplitBed(cython.str bedpath):
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


@cython.returns(cython.str)
def MergeBamList(bamlist, picardPath="default", memStr="-Xmx6G",
                 outbam="default"):
    """
    Merges a list of BAMs. Used for merging discordant read bams for
    parallelized VCF calls.
    """
    if(isinstance(bamlist, str)):
        bamlist = bamlist.split(":")
    if(outbam == "default"):
        outbam = bamlist[0].split(".")[0:-1] + ".merged.bam"
    commandStr = "java %s -jar %s AS=true" % (memStr, picardPath)
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
            raise ThisIsMadness(
                "Resubmission limit %s reached!" % self.maxresubs)
        self.popen = DevNullPopen(self.commandString)
        self.resubmissions += 1


@cython.returns(cython.str)
def GetOutVCFFromBMFsnvPopen(d):
    return d.commandString.split(" ")[9]


@cython.returns(cython.str)
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
    def __init__(self, stringlist, threads=4, MaxResubmissions=50,
                 func=None, cleanup=None):
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
            raise ThisIsMadness("func must be callable!")
        self.getReturnValue = func
        self.bgStrs = []
        self.fgStrs = []
        if(cleanup is not None):
            self.cleanup = cleanup

    @cython.returns(cython.long)
    def getJobNumber(self):
        return self.submitted + 1

    def submit(self):
        """
        Function for submitting bg (background) jobs.
        """
        if(len(self.queue) == 0):
            print("All jobs already submitted. :)")
            return
        while(len(self.dispatches) < self.threadcount):
            if(len(self.queue) != 0):
                print("Submitting job number %s." % str(self.getJobNumber()))
                self.submitted += 1
                cStr = self.queue.popleft()
                self.bgStrs.append(cStr)
                submitted = PopenCall(cStr)
                print("Command String: %s" % submitted.commandString)
                self.dispatches.append(submitted)
            else:
                print("All jobs submitted - check in later.")
                pass

    def check(self):
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
        self.submit()
        while len(self.queue) != 0:
            print("Submitting set of jobs for daemon.")
            self.submit()
            newqueue = tee(self.queue, 1)[0]
            try:
                nextStr = newqueue.next()
            except StopIteration:
                time.sleep(5)
                continue
            if("#" in nextStr):
                fgCStr = self.queue.popleft()
                self.fgStrs.append(fgCStr)
                print("Foreground submitting job #%s" % self.getJobNumber())
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
                time.sleep(5)
            self.check()
        print("All jobs submitted! Yay.")
        while(self.check() > 0 and len(self.queue) != 0):
            time.sleep(5)
            threadcount = self.check()
            if(threadcount < self.threadcount and len(self.queue) != 0):
                self.submit()
        while(self.check() != 0):
            """
            while(sum([d.popen.poll() is not None for d in
                       self.dispatches]) < len(self.dispatches)):
                pl("Sleeping because some jobs haven't returned yet.")
                time.sleep(5)
            """
            pl("Sleeping because some jobs haven't returned yet.")
            time.sleep(5)
        for key in self.outstrs.iterkeys():
            if(self.outstrs[key] is None):
                print("fgStrs: %s" % ":".join(self.fgStrs))
                print("bgStrs: %s" % ":".join(self.bgStrs))
                print(repr(self.outstrs))
                raise ThisIsMadness(
                    "Command string %s has a key value of None!" % key)
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
        if(not os.path.isfile(bam + ".bai")):
            raise ThisIsMadness("Index for bam not created!")
    return ziplist


def SplitBamByBedPysam(bampath, bedpath):
    """
    Rather than standard split by bed which uses parallel shell calls,
    instead this does it manually via pysam. Not sure if it's faster or not.
    """
    cdef pysam.calignmentfile.AlignedSegment rec
    cdef pysam.calignmentfile.AlignmentFile inHandle
    cdef cython.str bam
    pl("Getting bamlist.")
    bamlist = map(oig1, GetBamBedList(bampath, bedpath))
    refContigNumList = [ChrToPysamDict[bam.split(".")[-2]] for bam in bamlist]
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
        if(os.path.isfile(bam + ".bai") is False):
            raise ThisIsMadness("samtools couldn't index this bam. "
                                "It should already be sorted. Abort!")
    return bamlist


@cython.returns(cython.str)
def BMFsnvCommandString(cython.str bampath, cython.str conf="default",
                        cython.str bed="default"):
    """
    Returns a command string for calling bmftools snv
    """
    if(conf == "default"):
        raise ThisIsMadness("conf must be set for BMFsnvCommandString.")
    tag = str(uuid.uuid4().get_hex().upper()[0:8])
    if(bed == "default"):
        raise ThisIsMadness("bed must be set for BMFsnvCommandString.")
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


@cython.returns(cython.str)
def GetUUIDFromCommandString(cython.str cStr):
    """
    Gets the UUID tag from slave jobs. Used for concatenating VCFs from
    parallel calls.
    """
    return cStr.split("#G~")[1]


def GetBMFsnvPopen(bampath, bedpath, conf="default", threads=4,
                   parallel=False):
    """
    Makes a PopenDispatcher object for calling these variant callers.
    """
    if conf == "default":
        raise ThisIsMadness("conf file but be set for GetBMFsnvPopen")
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


@cython.returns(cython.str)
def TrimExt(cython.str fname):
    """
    Trims the extension from the filename so that I don't have to type this
    every single time.
    """
    return ".".join(fname.split(".")[:-1])


@cython.returns(cython.str)
def NPadSequence(cython.str seq, cython.long n=300):
    """
    Pads a sequence with "n" Ns.
    """
    return "N" * n + seq + "N" * n


@cython.returns(cython.str)
def FastaStyleSequence(cython.str seq):
    return ">" + seq + "\n" + seq


@cython.returns(cython.str)
def PadAndMakeFasta(cython.str seq, cython.long n=300):
    return FastaStyleSequence(NPadSequence(seq, n=n))


@cython.returns(cython.str)
def SequenceToFakeFq(cython.str seq):
    return ("@" + seq + "\n" + seq +
            "\n+\n" + "G" * len(seq))


@cython.returns(list)
def GetKmersToCheck(cython.str ref, cython.long k=30, list bedline=[],
                    cython.long padding=-1):
    """
    Gets a list of kmers which provide unique mappability
    to the region of interest.
    bedline should be the first 3 columns from a line in a bed file.
    """
    cdef cython.long i, start, end
    cdef list kmerList, candidateKmers
    if(padding < 0):
        pl("Padding not set - defaults to kmer size.")
        padding = k
    kmerList = []
    refHandle = pysam.FastaFile(ref)
    contig, start = bedline[0], bedline[1] - padding
    end = bedline[2] + padding
    regionStr = refHandle.fetch(contig, start, end)
    return [regionStr[i:i + k] for i in xrange(end - start - k)]


@cython.returns(cython.str)
def FastqStrFromKmerList(list kmerList):
    """
    Creates a dummy fastq string from a list of kmers.
    """
    return "\n".join(map(SequenceToFakeFq, kmerList))


@cython.returns(cython.str)
def BowtieFqToStr(cython.str fqStr, cython.str ref=None,
                  cython.long mismatches=-1, cython.long seed=-1):
    """
    Returns the string output of a bowtie call.
    """
    if(mismatches < 0):
        raise ThisIsMadness("mismatches must be set for BowtieFqToStr.")
    if(seed < 0):
        raise ThisIsMadness("seed length must be set for BowtieFqToStr.")
    cStr = ("echo %s | bowtie %s -a -n %s -l %s %s -" % (ref, mismatches,
                                                         seed, fqStr))
    return check_output(cStr, shell=True)


@cython.returns(list)
def GetMQPassReads(cython.str bwtStr, cython.long minMQ=1):
    """
    Takes a string output from bowtie and gets the names of the reads
    with MQ >= minMQ. Defaults to 1 (for a unique alignment)
    """
    cdef list readnames
    cdef cython.str line
    cdef tuple nameCount
    readnames = [line.split("\t")[0] for line in bwtStr.split("\n") if
                 line[0] != "@" and int(line.split("\t")[4]) >= minMQ]
    uniquereadnames = [nameCount[0] for nameCount in
                       cyfreq(readnames).iteritems() if
                       nameCount[1] == 1]
    pl("# of passing read names: %s" % len(uniquereadnames))
    return uniquereadnames


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

    def __init__(self, cython.str contig, cython.long start=-666,
                 cython.long end=-1, cython.long type=-137,
                 cython.str seq=None):
        self.contig = contig
        self.start = start
        self.end = end
        self.type = type
        self.seq = seq
        self.readnames = []
        self.uniqStr = None

    def __str__(self):
        raise ThisIsMadness("Abstract method must be inherited. Sorry, cdef w"
                            "on't let me actually make this an abstract "
                            "class.")

    @cython.returns(cython.long)
    def __len__(self):
        """
        Returns the number of reads supporting it which have been queried
        against this object since its creation.
        """
        return len(self.readnames)

    @cython.returns(cython.bint)
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

    @cython.returns(cython.str)
    def __getitem__(self, cython.long index):
        return self.readnames[index]

    @cython.returns(list)
    def sort(self):
        self.readnames = sorted(self.readnames)


cdef class Insertion(AbstractIndelContainer):
    """
    Expansion of the AbstractIndelContainer object.
    "Start" is actually the preceding base to the insertion (counting up),
    "End" is the following. end should always be greater than start.
    The important thing is to be able to hold read names and have a unique
    string representing each indel so that we can make calls.

    If no handle is provided, shen (Shannon Entropy) is set to -1.
    """

    def __init__(self, cython.str contig, cython.long start=-1,
                 cython.str seq=None, cython.str readname=None,
                 pysam.cfaidx.FastaFile handle=None,
                 cython.long window=20):
        if(start < 0):
            raise ThisIsMadness("start required for InsertionContainer.")
        self.contig = contig
        self.start = start
        self.end = start + 1
        self.type = 1
        self.seq = seq
        self.readnames = [readname]
        self.uniqStr = "Insertion|%s:%s,%s|%s" % (self.contig, self.start,
                                                  self.end, self.seq)
        try:
            self.shen = min([shen(handle.fetch(read.reference_id,
                                               start=start - window,
                                               end=start)),
                             shen(handle.fetch(read.reference_id,
                                               start=end, end=end + window))])
        except AttributeError:
            self.shen = -1.
        self.handle = handle
        self.shenwindow = window

    @cython.returns(cython.str)
    def __str__(self):
        return self.uniqStr + "|%s" % len(self.readnames)



cdef class Deletion(AbstractIndelContainer):
    """
    Expansion of the AbstractIndelContainer object.
    "Start" is the first missing base.
    "End" is the last missing base.
    If the deletion is of length one, start and end should be the same.
    """

    def __init__(self, cython.str contig, cython.long start=-1,
                 cython.long end=-1, cython.str readname=None,
                 pysam.cfaidx.FastaFile handle=None, cython.long window=20):
        self.contig = contig
        self.start = start
        self.end = end
        self.type = -1
        self.readnames = [readname]
        self.uniqStr = "Deletion|%s:%s,%s" % (self.contig, self.start,
                                              self.end)
        try:
            self.shen = min([shen(handle.fetch(read.reference_id,
                                               start=start - window,
                                               end=start)),
                             shen(handle.fetch(read.reference_id,
                                               start=end, end=end + window))])
        except AttributeError:
            self.shen = -1.
        self.handle = handle
        self.shenwindow = window

    @cython.returns(cython.str)
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
    def __init__(self, cython.str ref=None, cython.long window=10):
        self.data = {}
        self.counts = {}
        self.deletions = []
        self.insertions = []
        self.complexindels = []
        self.fastaRef = pysam.FastaFile(ref)
        self.window = window

    @cython.returns(cython.long)
    def __len__(self):
        return len(self.data)

    @cython.returns(list)
    def __getitem__(self, cython.str key):
        return self.data[key]

    def __setitem__(self, cython.str key, list value):
        self.data[key] = value

    def setIndelShen(self, AbstractIndelContainer_t indelObj):
        indelObj.shen = min([shen(self.fastaRef(indelObj.contig,
                                                indelObj.start - self.window,
                                                indelObj.start)),
                             shen(self.fastaRef(indelObj.contig, indelObj.end,
                                                indelObj.end + self.window))])
        indelObj.shenwindow = window

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

    def addIndel(self, AbstractIndelContainer_t indelObj,
                 cython.bint newWindowSize=False):
        try:
            self[indelObj.uniqStr] += indelObj.readnames
        except KeyError:
            self[indelObj.uniqStr] = indelObj.readnames
        self.setIndelShen(indelObj)
        if(isinstance(indelObj, Deletion)):
            self.deletions.append(indelObj)
        elif(isinstance(indelObj, Insertion)):
            self.insertions.append(indelObj)
        else:
            self.complexindels.append(indelObj)
        self.counts = {key: len(values) for key, values in
                       self.data.itervalues()}


    def mergeQuiver(self, IndelQuiver_t quiverObj):
        for key in quiverObj:
            try:
                self[key] += quiverObj[key]
            except KeyError:
                self[key] = quiverObj[key]
        self.counts = {key: len(values) for key, values in
                       self.data.itervalues()}
        self.deletions += quiverObj.deletions
        self.insertions += quiverObj.insertions
        self.complexindels += quiverObj.complexindels
        if(newWindowSize):
            self.window = quiverObj.window
        # Recalculate the Shannon entropy for each indel with new window
        for indel in self.deletions + self.insertions + self.complexindels:
            self.setIndelShen(indel)

    @cython.returns(dict)
    def getIndelCounter(self, cython.str uniqStr):
        try:
            return cyfreq(self[uniqStr])
        except KeyError:
            return {}


def hamming_cousins_exact(cython.str s, cython.long n,
                          set alphabet={"A", "C", "G", "T"}):
    """Generate strings over alphabet whose Hamming distance from s is
    exactly n.

    >>> sorted(hamming_cousins_exact('abc', 0, 'abc'))
    ['abc']
    >>> sorted(hamming_cousins_exact('abc', 1, 'abc'))
    ['aac', 'aba', 'abb', 'acc', 'bbc', 'cbc']
    >>> sorted(hamming_cousins_exact('aaa', 2, 'ab'))
    ['abb', 'bab', 'bba']

    """
    cdef cython.long lalphabetSub1
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


def hamming_cousins(cython.str s, cython.long n,
                    set alphabet={"A", "C", "G", "T"}):
    """Generate strings over alphabet whose Hamming distance from s is
    less than or equal to n.

    >>> sorted(hamming_cousins('abc', 0, 'abc'))
    ['abc']
    >>> sorted(hamming_cousins('abc', 1, 'abc'))
    ['aac', 'aba', 'abb', 'abc', 'acc', 'bbc', 'cbc']
    >>> sorted(hamming_cousins('aaa', 2, 'ab'))
    ['aaa', 'aab', 'aba', 'abb', 'baa', 'bab', 'bba']

    """
    return chain(*(hamming_cousins_exact(s, i, alphabet) for i in range(n + 1)))


@cython.returns(Insertion_t)
def GetInsertionFromAlignedSegment(pysam.calignmentfile.AlignedSegment read,
                                   pysam.cfaidx.FastaFile handle=None):
    """
    Creates an Insertion object under the assumption that there is
    only one set of inserted bases in the read. I should come back and
    generalize it...
    """
    cdef tuple InsInf
    InsInf = GetInsertedStrs(read)[0]
    return Insertion(PysamToChrDict[read.reference_id], start=InsInf[1],
                     seq=InsInf[0], readname=read.query_name,
                     handle=handle)


@cython.returns(Deletion_t)
def GetDeletionFromAlignedSegment(pysam.calignmentfile.AlignedSegment read,
                                  pysam.cfaidx.FastaFile handle=None):
    """
    Creates a Deletion object from a read.
    """
    cdef cython.long start, end
    cdef list coords
    coords = GetDeletedCoordinates(read)
    start = coords[0]
    end = coords[-1]
    return Deletion(PysamToChrDict[read.reference_id], start=start,
                    end=end, readname=read.query_name, handle=handle)
