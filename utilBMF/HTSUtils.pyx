# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True
import logging
import shlex
import subprocess
from collections import Counter
from itertools import groupby
from operator import itemgetter
from operator import attrgetter
import copy
import uuid
import os
import operator

import pysam
import numpy as np
cimport numpy as np
from Bio.Seq import Seq
import cython
import numconv

import MawCluster
from utilBMF.ErrorHandling import *


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
    print(string.replace(
          "\t", "\\t").replace("\n", "\\n").replace(
          "'", "\'").replace('"', '\\"'))
    return


def ToStr(x):
    """
    Needed to use map for str conversion
    """
    return str(x)


# TODO: Write something to create these dictionaries from a SAM header


PysamToChrDict = {}
for i in range(22):
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
for i in range(22):
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

"""
This dictionary is required by some tools in the package to deal with the fact
that pysam doesn't store contig names, only which contig number it is,
based on the order found in the sam header.
"""


DefaultSamHeader = repr("{'SQ': [{'LN': 249250621, 'SN': '1'}, {'LN': 2431"
                        "99373, 'SN': '2'}, {'LN': 198022430, 'SN': '3'}, "
                        "{'LN': 191154276, 'SN': '4'}, {'LN': 180915260, '"
                        "SN': '5'}, {'LN': 171115067, 'SN': '6'}, {'LN': 1"
                        "59138663, 'SN': '7'}, {'LN': 146364022, 'SN': '8'"
                        "}, {'LN': 141213431, 'SN': '9'}, {'LN': 135534747"
                        ", 'SN': '10'}, {'LN': 135006516, 'SN': '11'}, {'L"
                        "N': 133851895, 'SN': '12'}, {'LN': 115169878, 'SN"
                        "': '13'}, {'LN': 107349540, 'SN': '14'}, {'LN': 1"
                        "02531392, 'SN': '15'}, {'LN': 90354753, 'SN': '16"
                        "'}, {'LN': 81195210, 'SN': '17'}, {'LN': 78077248"
                        ", 'SN': '18'}, {'LN': 59128983, 'SN': '19'}, {'LN"
                        "': 63025520, 'SN': '20'}, {'LN': 48129895, 'SN': "
                        "'21'}, {'LN': 51304566, 'SN': '22'}, {'LN': 15527"
                        "0560, 'SN': 'X'}, {'LN': 59373566, 'SN': 'Y'}, {'"
                        "LN': 16569, 'SN': 'MT'}, {'LN': 4262, 'SN': 'GL00"
                        "0207.1'}, {'LN': 15008, 'SN': 'GL000226.1'}, {'LN"
                        "': 19913, 'SN': 'GL000229.1'}, {'LN': 27386, 'SN'"
                        ": 'GL000231.1'}, {'LN': 27682, 'SN': 'GL000210.1'"
                        "}, {'LN': 33824, 'SN': 'GL000239.1'}, {'LN': 3447"
                        "4, 'SN': 'GL000235.1'}, {'LN': 36148, 'SN': 'GL00"
                        "0201.1'}, {'LN': 36422, 'SN': 'GL000247.1'}, {'LN"
                        "': 36651, 'SN': 'GL000245.1'}, {'LN': 37175, 'SN'"
                        ": 'GL000197.1'}, {'LN': 37498, 'SN': 'GL000203.1'"
                        "}, {'LN': 38154, 'SN': 'GL000246.1'}, {'LN': 3850"
                        "2, 'SN': 'GL000249.1'}, {'LN': 38914, 'SN': 'GL00"
                        "0196.1'}, {'LN': 39786, 'SN': 'GL000248.1'}, {'LN"
                        "': 39929, 'SN': 'GL000244.1'}, {'LN': 39939, 'SN'"
                        ": 'GL000238.1'}, {'LN': 40103, 'SN': 'GL000202.1'"
                        "}, {'LN': 40531, 'SN': 'GL000234.1'}, {'LN': 4065"
                        "2, 'SN': 'GL000232.1'}, {'LN': 41001, 'SN': 'GL00"
                        "0206.1'}, {'LN': 41933, 'SN': 'GL000240.1'}, {'LN"
                        "': 41934, 'SN': 'GL000236.1'}, {'LN': 42152, 'SN'"
                        ": 'GL000241.1'}, {'LN': 43341, 'SN': 'GL000243.1'"
                        "}, {'LN': 43523, 'SN': 'GL000242.1'}, {'LN': 4369"
                        "1, 'SN': 'GL000230.1'}, {'LN': 45867, 'SN': 'GL00"
                        "0237.1'}, {'LN': 45941, 'SN': 'GL000233.1'}, {'LN"
                        "': 81310, 'SN': 'GL000204.1'}, {'LN': 90085, 'SN'"
                        ": 'GL000198.1'}, {'LN': 92689, 'SN': 'GL000208.1'"
                        "}, {'LN': 106433, 'SN': 'GL000191.1'}, {'LN': 128"
                        "374, 'SN': 'GL000227.1'}, {'LN': 129120, 'SN': 'G"
                        "L000228.1'}, {'LN': 137718, 'SN': 'GL000214.1'}, "
                        "{'LN': 155397, 'SN': 'GL000221.1'}, {'LN': 159169"
                        ", 'SN': 'GL000209.1'}, {'LN': 161147, 'SN': 'GL00"
                        "0218.1'}, {'LN': 161802, 'SN': 'GL000220.1'}, {'L"
                        "N': 164239, 'SN': 'GL000213.1'}, {'LN': 166566, '"
                        "SN': 'GL000211.1'}, {'LN': 169874, 'SN': 'GL00019"
                        "9.1'}, {'LN': 172149, 'SN': 'GL000217.1'}, {'LN':"
                        " 172294, 'SN': 'GL000216.1'}, {'LN': 172545, 'SN'"
                        ": 'GL000215.1'}, {'LN': 174588, 'SN': 'GL000205.1"
                        "'}, {'LN': 179198, 'SN': 'GL000219.1'}, {'LN': 17"
                        "9693, 'SN': 'GL000224.1'}, {'LN': 180455, 'SN': '"
                        "GL000223.1'}, {'LN': 182896, 'SN': 'GL000195.1'},"
                        " {'LN': 186858, 'SN': 'GL000212.1'}, {'LN': 18686"
                        "1, 'SN': 'GL000222.1'}, {'LN': 187035, 'SN': 'GL0"
                        "00200.1'}, {'LN': 189789, 'SN': 'GL000193.1'}, {'"
                        "LN': 191469, 'SN': 'GL000194.1'}, {'LN': 211173, "
                        "'SN': 'GL000225.1'}, {'LN': 547496, 'SN': 'GL0001"
                        "92.1'}], 'PG': [{'PN': 'bwa', 'ID': 'bwa', 'VN': "
                        "'0.7.10-r789', 'CL': 'bwa mem -t 4 -v 1 -Y -M -T "
                        "0'}], 'HD': {'SO': 'queryname', 'VN': '1.4'}}")


def GetPysamToChrDict(alignmentFileText):
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
                      PysamToChrDict.keys()}
    return


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
    raise ThisIsMadness(string)


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
        minLen = readPair.read1.query_length / 2
    else:
        raise ThisIsMadness("minShare parameter required. Integer for "
                            "an absolute number of bases overlapped re"
                            "quired, float for a fraction of read length.")
    return sum([x == 2 for x in
                Counter(readPair.read1.get_reference_positions() +
                        readPair.read2.get_reference_positions()).values()]
               ) >= minLen


def align_bwa_aln(R1, R2, ref="default", opts="", outBAM="default"):
    """
    Aligns a set of paired-end reads using bwa aln. Defaults to 4 threads.
    """
    if(ref == "default"):
        FacePalm("Reference file index required for alignment!")
    if(opts == ""):
        opts = "-n 3 -t 4"
    if(outBAM == "default"):
        outBAM = '.'.join(R1.split('.')[0:-1]) + ".aln.bam"
    str(uuid.uuid4().get_hex().upper()[0:8])
    R1Sai = R1 + ".tmp.sai"
    R2Sai = R2 + ".tmp.sai"
    alnStr1 = ("bwa aln " + opts + " " + " ".join([ref, R1]) +
               " > " + R1Sai)
    PipedShellCall(alnStr1)
    alnStr2 = ("bwa aln " + opts + " " + " ".join([ref, R2]) +
               " > " + R2Sai)
    PipedShellCall(alnStr2)
    sampeStr = ("bwa sampe " " " +
                " ".join([ref, R1Sai, R2Sai, R1, R2]) +
                " | samtools view -h - > " + outBAM)
    printlog("bwa aln string: {}".format(sampeStr))
    PipedShellCall(sampeStr)
    os.remove(R1Sai)
    os.remove(R2Sai)
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
        command_str = ('bwa mem {} {} {} {}'.format(opt_concat, ref, R1, R2) +
                       " | samtools view -Sbh - > {}".format(outBAM))
    else:
        command_str = (path + ' mem {} {} {}'.format(opt_concat, ref, R1, R2) +
                       " {} | samtools view -Sbh - > {}".format(R2, outBAM))
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
    if(delete is True):
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
    # if(isinstance(bedRef, str) is True):
    #     bedRef = ParseBed(bedRef)
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
    if(isinstance(bedRef, str) is True):
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
                        uuid="true",
                        threads="4"):
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
                  " -@ {} {}".format(threads, inBAM))
    printlog("About to call sort command: {}".format(CommandStr))
    subprocess.check_call(shlex.split(CommandStr))
    printlog("Now indexing.")
    subprocess.check_call(shlex.split("samtools index {}".format(outBAM)))
    return outBAM


def NameSort(inBAM, outBAM="default", prefix="MetasyntacticVar",
             uuid="true", threads="4"):
    # If uuid is either a boolean true or is a string containing true,
    # then a random string is generated for the output
    if(str(uuid).lower() == "true"):
        import uuid
        prefix += str(uuid.uuid4().get_hex().upper()[0:8])
    if(outBAM == "default"):
        outBAM = '.'.join(inBAM.split('.')[0:-1]) + '.NameSort.bam'
    CommandStr = ("samtools sort -T {} -O bam -o {}".format(prefix, outBAM) +
                  " -@ {} -n {}".format(threads, inBAM))
    printlog("About to call sort command: {}".format(CommandStr))
    subprocess.check_call(shlex.split(CommandStr))
    printlog("Namesort successful, sorted bam available at: {}".format(outBAM))
    return outBAM


def mergeBam(samList, memoryStr="-XmX16",
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


class ReadPair:

    """
    Holds both bam record objects in a pair.
    Currently, one read unmapped and one read soft-clipped are
    both marked as soft-clipped reads.
    """

    def __init__(self, read1, read2):
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
            if("S" in read1.cigarstring is True):
                self.read1_soft_clipped = True
            else:
                self.read1_soft_clipped = False
        if(read2.is_unmapped):
            self.read2_is_unmapped = True
            self.read2_soft_clipped = True
        else:
            self.read2_is_unmapped = False
            if("S" in read2.cigarstring is True):
                self.read2_soft_clipped = True
            else:
                self.read2_soft_clipped = False
        if(self.read1_soft_clipped is True):
            self.read1_softclip_seqs = []
        if(self.read2_soft_clipped is True):
            self.read2_softclip_seqs = []
        if(self.read1_is_unmapped is True):
            self.read1_contig = "*"
        else:
            self.read1_contig = PysamToChrDict[read1.reference_id]
        if(self.read2_is_unmapped is True):
            self.read2_contig = "*"
        else:
            self.read2_contig = PysamToChrDict[read2.reference_id]
        self.SameContig = (read1.reference_id == read2.reference_id)
        self.ContigString = ",".join(sorted([self.read1_contig,
                                             self.read2_contig]))
        # TODO: write a script to create an array of soft-clipped sequences
        # from each read

    def NumOverlappingBed(self, bedLines="default"):
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

    def getReads(self):
        return [self.read1, self.read2]


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
    if(np.any(np.less(Pair.read1.query_qualities, minBQ)) or
       np.any(np.less(Pair.read2.query_qualities, minBQ))):
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
        return sorted(RecordsArray, key=abs(attrgetter("insert_size")))
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


def ParseBed(bedfile):
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


def parseConfig(string):
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
    ConfigDict = {}
    for line in parsedLines:
        ConfigDict[line.split("=")[0].strip()] = line.split("=")[1].strip()
    return ConfigDict


def ReadListToCovCounter(reads, minClustDepth=3, minPileupLen=10):
    """
    Makes a Counter object of positions covered by a set of reads.
    Only safe at this point for intrachromosomal rearrangements!
    """
    return Counter(reduce(lambda x, y: x + y,
                          [r.get_reference_positions() for r in reads]))


def ReadPairListToCovCounter(ReadPairList, minClustDepth=5, minPileupLen=10):
    """
    Makes a Counter object of positions covered by a set of read pairs.
    Only safe at this point for intrachromosomal rearrangements!
    We discount the "duplex" positions because we want to look for pileups of
    read pairs (ultimately, for supporting a structural variant).
    """
    posList = []
    posListDuplex = []
    for pair in ReadPairList:
        R1Pos = pair.read1.get_reference_positions()
        R2Pos = pair.read2.get_reference_positions()
        operator.iadd(posList, R1Pos)
        operator.iadd(posList, R2Pos)
        operator.iadd(posListDuplex, [pos for pos in R1Pos if pos in R2Pos])
    PosCounts = Counter(posList)
    PosDuplexCounts = Counter(posListDuplex)
    # decrement the counts for each position to account for
    # both reads in a pair mapping to the same location.
    for key in PosDuplexCounts.keys():
        PosCounts[key] -= PosDuplexCounts[key]
    PosCounts = dict([i for i in PosCounts.items()
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


def SplitSCRead(read):
    try:
        assert isinstance(read, pysam.calignmentfile.AlignedSegment)
    except AssertionError:
        FacePalm("You can't split a read by soft-clipping if it's not a read!")
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


def CreateIntervalsFromCounter(CounterObj, minPileupLen=0, contig="default",
                               bedIntervals="default",
                               mergeDist=0, minClustDepth=5):
    """
    From a Counter object containing the sum of the output of
    get_reference_positions for a list of AlignedSegment objects, it creates a
    list of continuous intervals, calculates the mean coverage for each, and
    returns two lists, the first containing the 0-based open intervals and
    second containing the mean coverage of that interval.
    bedIntervals must be in ParseBed output format.
    """
    assert isinstance(CounterObj, dict)
    CounterObj = dict([i for i in CounterObj.items()
                      if i[1] >= minClustDepth])
    IntervalList = []
    MeanCovList = []
    if(contig == "default"):
        FacePalm("contig required for this function!")
    for k, g in groupby(
            enumerate(sorted(CounterObj.keys())),
            lambda i_x: i_x[0] - i_x[1]):
        posList = map(itemgetter(1), g)
        if(posList[0] < posList[-1]):
            interval = [contig, posList[0], posList[-1] + 1]
        else:
            interval = [contig, posList[-1], posList[0] + 1]
        if(interval[2] - interval[1] < minPileupLen):
            print("Interval too short. Length: {}".format(interval[2]
                                                          - interval[1]))
            continue
        IntervalList.append(interval)
        MeanCovList.append(np.mean([CounterObj[key] for key in posList if
                                    len(posList) >= minPileupLen]))
    # Now merge, in case some are entirely adjacent
    if(len(IntervalList) == 0):
        return []
    print("Now attempting to merge any adjacent intervals. Number: {}".format(
        len(IntervalList)))
    IntervalList = sorted(IntervalList, key=itemgetter(1))
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
for i in range(100):
    ph2chrDict[i] = chr(i + 33)
chr2ph = {ph2chrDict[key]: key for key in
          ph2chrDict.keys()}

"""
@cython.returns(np.int64_t)
def chr2ph(x):
    \"""
    Converts a character to its corresponding phred integer representation
    \"""
    return ord(x) - 33
"""


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
    c = map(operator.itemgetter(1), cigar)
    cumSum = [sum(c[:i + 1]) for i in range(len(c))]
    for n, entry in enumerate(cigar):
        if n == 0:
            tuples.append((entry[0], range(entry[1])))
        else:
            tuples.append((entry[0], range(cumSum[n - 1], cumSum[n])))
    return tuples


@cython.locals(cigarOp=cython.int)
def GetQueryIndexForCigarOperation(read, cigarOp=-1):
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


def GetReadNucsForFiltCigar(filtCigar, read):
    """
    Returns a list of strings for nucleotides matching a filtCigar.
    """
    seq = read.seq
    return ["".join([seq[i] for i in g[1]]) for g in filtCigar]


def GetGenomicCoordsForFiltCigar(filtCigar, read):
    """
    Returns a list of lists for genomic coordinates matching a filtCigar.
    """
    dgap = dict(read.get_aligned_pairs())
    return [[dgap[i] for i in g[1]] for g in filtCigar]


def GetReadSequenceForCigarOp(read, cigarOp=-1):
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
    return


def GetDeletedCoordinates(read):
    """
    Returns a list of integers of genomic coordinates for deleted bases.
    """
    assert isinstance(read, pysam.calignmentfile.AlignedSegment)

    return sorted([i[1] for i in read.get_aligned_pairs() if i[0] is None])
