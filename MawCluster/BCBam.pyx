#!python
# cython: c_string_type=str, c_string_encoding=ascii
# cython: boundscheck=False
import shlex
import subprocess
import os
import shutil
import logging
from copy import copy as ccopy
from os import path
import operator
from operator import attrgetter as oag, methodcaller as mc
import string
import uuid
import sys
from subprocess import check_call
from functools import partial
from itertools import groupby, chain
from array import array

import numpy as np
from numpy import (sum as nsum, multiply as nmul,
                   subtract as nsub, argmax as nargmax,
                   vstack as nvstack, char)
from cytoolz import frequencies as cyfreq
import pysam
import cython

from .BCFastq import GetDescriptionTagDict as getdesc
from . import BCFastq
from utilBMF.HTSUtils import (printlog as pl,
                              FractionAligned, FractionSoftClipped,
                              SWRealignAS, pPileupRead, BedtoolsBam2Fq,
                              BwaswCall, samtoolsMergeBam, pFastqProxy,
                              TrimExt)
from utilBMF.ErrorHandling import (IllegalArgumentError, ThisIsMadness as Tim,
                                   MissingExternalTool)
from utilBMF import HTSUtils
from warnings import warn
import SecC


def AbraCadabra(inBAM, outBAM="default",
                jar="default", sortMem="default", ref="default",
                threads="4", bed="default", working="default",
                log="default", bint fixMate=True, tempPrefix="tmpPref",
                rLen=-1, intelPath="default", bint leftAlign=True,
                bint kmers_precomputed=True):
    """
    Calls abra for indel realignment. It supposedly
    out-performs GATK's IndelRealigner, though it does right-align
    some indels.

    It also calls samtools fixmate to restore mate information and
    bamleftalign to left align any that abra right-aligned.
    """
    if(rLen < 0):
        raise IllegalArgumentError("Read length must be set to call abra due"
                                   " to the benefits of inferring ideal para"
                                   "meters from the !")
    if(jar == "default"):
        raise MissingExternalTool("Required: Path to abra jar!")
    else:
        pl("Non-default abra jar used: " + jar)
    if(sortMem == "default"):
        sortMem = "-Xmx16G"
        pl("Default memory string used: " + sortMem)
    else:
        pl("Non-default memory string used: " + sortMem)
    if(ref == "default"):
        raise ValueError("Reference fasta must be provided!")
    if(ref.split(".")[-1] == "gz"):
        warn("Reference fasta is gzipped, with which "
             "abra is not compatible. Be warned!", UserWarning)
    else:
        pl("Reference file set: {}.".format(ref))
    if(bed == "default"):
        raise ValueError("Bed file required.")
    else:
        pl("Bed file set: {}.".format(bed))
    if(working == "default"):
        bamFilename= path.basename(inBAM)
        working = (path.dirname(inBAM) + bamFilename.split('.')[0] +
                   ".working_dir")
        pl("Default working directory set to be: " + working)
    else:
        pl("Non-default working directory: " + working)
    if(log == "default"):
        log = "abra.log"
    if(outBAM == "default"):
        outBAM = '.'.join(inBAM.split('.')[0:-1]) + '.abra.bam'
    pl(("Command to reproduce the call of this function: "
        "AbraCadabra(\"{}\", outBAM=\"{}\", jar=\"{}\", ".format(inBAM,
                                                                 outBAM,
                                                                 jar) +
        "sortMem=\"{}\", ref=\"{}\", threads=\"{}\", ".format(sortMem,
                                                             ref, threads) +
        "bed=\"{}\", working=\"{}\", log=\"{}\")".format(bed, working, log)))
    if(path.isdir(working)):
        pl("Working directory already exists - deleting!")
        shutil.rmtree(working)
    # Check bed file to make sure it is in appropriate format for abra
    if(kmers_precomputed is False):
        bed = AbraKmerBedfile(bed, ref=ref, abra=jar,
                              rLen=rLen)
    if(path.isfile(inBAM + ".bai") is False):
        pl("No bam index found for input bam - attempting to create.")
        check_call(['samtools', 'index', inBAM])
        if(path.isfile(inBAM + ".bai") is False):
            inBAM = HTSUtils.CoorSortAndIndexBam(inBAM, outBAM, uuid=True)
    command = ("java {} -jar {} --in {}".format(sortMem, jar, inBAM) +
               " --out {} --ref {} --targets".format(outBAM, ref) +
               " {} --threads {} ".format(bed, threads) +
               "--working %s --mbq 200 --mer 0.0025 --mad 20000" % working)
    if(kmers_precomputed):
        command = command.replace("--targets", "--target-kmers")
    pl("Command: {}.".format(command))
    check_call(shlex.split(command), shell=False)
    pl("Deleting abra's intermediate directory.")
    check_call(["rm", "-rf", working])
    if(fixMate):
        pl("Now fixing mates after abra's realignment.")
        tempFilename = tempPrefix + str(
            uuid.uuid4().get_hex()[0:8]) + ".working.tmp"
        nameSorted = HTSUtils.NameSort(outBAM)
        commandStrFM = "samtools fixmate %s %s -O bam" % (nameSorted,
                                                          tempFilename)
        check_call(shlex.split(commandStrFM))
        check_call(["rm", "-rf", nameSorted])
        check_call(["mv", tempFilename, outBAM])
    if(leftAlign):
        # Calls bamleft align to make sure things are fixed up.
        tmpfile = str(uuid.uuid4().get_hex()[0:8]) + '.bam'
        cStr = ("samtools view -ubh %s | bamleftalign -f " % (outBAM) +
                "%s -c > %s && mv %s %s" % (ref, tmpfile, tmpfile, outBAM))
        check_call(cStr, shell=True)
    return outBAM


@cython.locals(rLen=int)
def AbraKmerBedfile(inbed, rLen=-1, ref="default", outbed="default",
                    nt=4, abra="default"):
    if(abra == "default"):
        raise MissingExternalTool(
            "Path to abra jar required for running KmerSizeCalculator.")
    if(ref == "default"):
        raise Tim(
            "Path to reference required for running KmerSizeCalculator.")
    if(inbed == "default"):
        raise Tim(
            "Path to input bed required for running KmerSizeCalculator.")
    if(rLen < 0):
        raise Tim(
            "Read length required for running KmerSizeCalculator.")
    if(outbed == "default"):
        outbed = ".".join(inbed.split(".")[0:-1] + ["abra", "kmer"])
    commandStr = ("java -cp %s abra.KmerSizeEvaluator " % abra +
                  "%s %s %s %s %s" % (rLen, ref, outbed, nt, inbed))
    pl("AbraKmerSizeEvaluator call string: %s" % commandStr)
    check_call(commandStr, shell=True)
    return outbed


def Bam2Sam(inBAM, outsam):
    pl("Bam2Sam. Input: {}. Output: {}.".format(inBAM, outsam))
    output = open(outsam, 'w', 0)
    command_str = 'samtools view -h {}'.format(inBAM)
    pl(command_str)
    check_call(shlex.split(command_str), stdout=output, shell=False)
    return(command_str, outsam)


def GATKIndelRealignment(inBAM, gatk="default", ref="default",
                         bed="default", dbsnp="default"):
    if(ref == "default"):
        raise MissingExternalTool("Reference file required"
                                  " for Indel Realignment")
    if(bed == "default"):
        raise Tim("Bed file required for Indel Realignment")
    if(gatk == "default"):
        raise MissingExternalTool("Path to GATK Jar required "
                                  "for Indel Realignment")
    print dbsnp
    if(dbsnp == "default"):
        dbsnpStr = ""
        pl("Running GATK Indel Realignment without dbSNP for known indels.")
    else:
        dbsnpStr = " -known %s " % dbsnp
    out = ".".join(inBAM.split(".")[0:-1] + ["realignment", "intervals"])
    outBAM = ".".join(inBAM.split(".")[0:-1] + ["gatkIndelRealign", "bam"])
    RTCString = "".join([
        "java -jar %s -T RealignerTargetCreator" % gatk,
        " -R %s -o %s -I %s -L:intervals,BED %s" % (ref, out, inBAM, bed),
        dbsnpStr])
    pl("RealignerTargetCreator string: %s" % RTCString)
    try:
        check_call(shlex.split(RTCString))
    except subprocess.CalledProcessError:
        pl("GATK RealignerTargetCreator failed. Still finish the "
           "analysis pipeline...")
        return inBAM
    IRCString = "".join(["java -jar %s -T IndelRealigner -targetInt" % gatk,
                         "ervals %s -R %s -I %s -o %s " % (out, ref,
                                                           inBAM, outBAM),
                         dbsnpStr])
    pl("IndelRealignerCall string: %s" % IRCString)
    try:
        check_call(shlex.split(IRCString))
    except subprocess.CalledProcessError:
        pl("GATK IndelRealignment failed. Still finish the analysis pipeline.")
        return inBAM
    pl("Successful GATK indel realignment. Output: %s" % outBAM)
    return outBAM


cdef dict cGetCOTagDict(AlignedSegment_t read):
    cdef cystr s, cStr
    cStr = read.opt("CO")
    return dict([s.split("=") for s in cStr.split("|")[1:]])


cpdef dict pGetCOTagDict(AlignedSegment_t read):
    return cGetCOTagDict(read)


@cython.wraparound(False)
cdef inline char * cRPString(bam1_t * src, bam_hdr_t * hdr) nogil:
    cdef char[100] buffer
    cdef size_t length
    length = sprintf("%s:%i,%s:%i", buffer, hdr.target_name[src.core.tid],
                     src.core.pos, hdr.target_name[src.core.mtid],
                     src.core.mpos)
    return buffer


cdef class BamPipe:
    """
    Creates a callable function which acts on a BAM stream.

    :param function - callable function which returns an input BAM object.
    :param bin_input - boolean - true if input is BAM
    false for TAM/SAM
    :param bin_output - boolean - true to output in BAM format.
    :param uncompressed_output - boolean - true to output uncompressed
    BAM records.
    """
    cpdef process(self):
        cdef AlignedSegment_t read
        [self.write(self.function(read)) for read in self.inHandle]

    def __init__(self, object function, bint bin_input, bint bin_output,
                 bint uncompressed_output=False):
        if(bin_input):
            self.inHandle = pysam.AlignmentFile("-", "rb")
        else:
            self.inHandle = pysam.AlignmentFile("-", "r")
        if(bin_output):
            if(uncompressed_output):
                self.outHandle = pysam.AlignmentFile(
                    "-", "wbu", template=self.inHandle)
            else:
                self.outHandle = pysam.AlignmentFile(
                    "-", "wb", template=self.inHandle)
        assert hasattr("__call__", function)
        self.function = function

    cdef write(self, AlignedSegment_t read):
        self.outHandle.write(read)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double getSF(AlignedSegment_t read):
    cdef tuple tup
    cdef int sum, sumSC
    sum = 0
    sumSC = 0
    if(read.cigarstring is None):
        return 0.
    for tup in read.cigar:
        sum += tup[1]
        if(tup[0] == 4):
            sumSC += tup[1]
    return sumSC * 1. / sum if(sum != 0) else 0.


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double getAF(AlignedSegment_t read):
    cdef tuple tup
    cdef int sum, sumAligned
    sum = 0
    sumAligned = 0
    if(read.cigarstring is None):
        return 0.
    for tup in read.cigar:
        sum += tup[1]
        if(tup[0] == 0):
            sumAligned += tup[1]
    return sumAligned * 1. / sum


cdef inline void CompareSeqQual(int32_t * template_quality,
                                int32_t * cmp_quality,
                                int8_t * seq1, int8_t * seq2,
                                int32_t rLen) nogil:
    cdef size_t index
    for index in range(rLen):
        if(seq1[index] == seq2[index]):
            template_quality[index] = MergeAgreedQualities(
                <int>template_quality[index], <int>cmp_quality[index]
            )
        else:
            if(template_quality[index] < cmp_quality[index]):
                seq1[index] = seq2[index]
            template_quality[index] = MergeDiscQualities(
                <int>template_quality[index], <int>cmp_quality[index]
            )
    return


cdef AlignedSegment_t MergeBamRecs(AlignedSegment_t template,
                                   AlignedSegment_t cmp):
    cdef py_array qual1, qual2, seq1, seq2
    qual1 = template.opt("PV")
    qual2 = cmp.opt("PV")
    seq1 = array('B')
    seq2 = array('B')
    c_array.resize(seq1, template._delegate.core.l_qseq)
    c_array.resize(seq2, template._delegate.core.l_qseq)
    memcpy(seq1.data.as_shorts, bam_get_seq(template._delegate),
           template._delegate.core.l_qseq)
    memcpy(seq2.data.as_shorts, bam_get_seq(cmp._delegate),
           template._delegate.core.l_qseq)
    CompareSeqQual(<int32_t *>qual1.data.as_ints,
                   <int32_t *>qual2.data.as_ints,
                   <int8_t *>seq1.data.as_shorts,
                   <int8_t *>seq2.data.as_shorts,
                   template._delegate.core.l_qseq)
    template.set_tag("PV", qual1)
    template.query_sequence = seq1.tostring()
    return template


cdef size_t build_dist_mtx(list recs, ndarray[int8_t, ndim=2] distmtx,
                           int8_t bLen):
    cdef size_t qcount, ccount, size, maxfm, FM, idxmaxfm
    cdef AlignedSegment_t qrec, cmprec
    size = len(recs)
    distmtx = np.zeros([size, size], dtype=np.int8)
    maxfm = 0
    for qcount, qrec in enumerate(recs):
        FM = qrec.opt("FM")
        if FM > maxfm:
            maxfm = FM
            idxmaxfm = qcount
        ccount = qcount + 1
        for cmprec in recs[ccount:]:
            distmtx[qcount,ccount] = pBarcodeHD(qrec, cmprec, bLen)
            ccount += 1
    return idxmaxfm


cdef list BFF(list recs, int8_t bLen, char mmlim):

    # C definitions
    cdef AlignedSegment_t qrec, cmprec
    cdef size_t size, qcount, ccount, idxmaxfm, FM, maxfm
    cdef list ret, merging_sets, sublist
    cdef ndarray[int8_t, ndim=2] distmtx
    cdef py_array mergeidx, fms
    cdef int8_t i

    # Debug message
    sys.stderr.write("Now attempting to flatten this set of"
                     " records which share coordinates.\n")

    # Setup
    size = len(recs)
    mergeidx = array('B')
    fms = array('i')
    c_array.resize(mergeidx, size)
    c_array.resize(fms, size)
    distmtx = np.zeros([size, size], dtype=np.int8)

    # Calculate the hamming distances. No way around that.
    idxmaxfm = 0
    maxfm = 0
    for qcount, qrec in enumerate(recs):
        # Note that I need to initialize these arrays, as I didn't set them
        # as clear when resizing.
        FM = qrec.opt("FM")
        fms[qcount] = FM
        if FM > maxfm:
            maxfm = FM
            idxmaxfm = qcount
        mergeidx[qcount] = 0
        ccount = qcount + 1

        for cmprec in recs[ccount:]:
            distmtx[qcount,ccount] = pBarcodeHD(qrec, cmprec, bLen)
            if distmtx[qcount,ccount] <= mmlim:
                mergeidx[qcount] = 1
            ccount += 1

    # Remove the ones which don't need to be flattened.
    ret = []
    qcount = size
    for qrec in recs[size:0:-1]:
        if not mergeidx[qcount]:
            ret.append(qrec)
            recs.remove(qrec)
        qcount -= 1

    # "Map" out which sets of reads need to be flattened into eah other.
    # First, get the one with the largest family size, then lump it into
    # a merging_sets sublist to be flattened.
    # and all within mmlim of it into a single record.
    merging_sets = []
    while True:
        size = len(recs)
        if size == 0:
            break
        idxmaxfm = build_dist_mtx(recs, distmtx, bLen)
        FilterSet(recs, distmtx, merging_sets, idxmaxfm, mmlim)

    ret += [MergeBamRecList(sublist) for sublist in merging_sets]
    sys.stderr.write("Now returning my list of records to write.\n")
    return ret


cdef void FilterSet(list recs,
                    ndarray[int8_t, ndim=1, mode="c"] distmtx,
                    list merging_sets, size_t idxmaxfm, char mmlim):
    # FilterSet needs to add the relevant reads to a new list in
    # merging_sets, and remove the relevant records from recs
    cdef size_t index
    cdef AlignedSegment_t rec
    cdef py_array indices
    cdef list to_append
    indices = array('i')
    to_append = []
    for index, rec in enumerate(recs):
        if idxmaxfm == index or distmtx[idxmaxfm,index] < mmlim:
            indices.append(index)
    for index in range(len(indices), 0, -1):
        rec = recs[index]
        to_append.append(rec)
        recs.remove(rec)
    merging_sets.append(to_append)
    return


cpdef inline cystr pBamRescue(cystr inBam, cystr outBam,
                              char mmlim, int8_t bLen):
    """
    :param inBam: [cystr/arg] - path to input bam
    :param outBam: [cystr/arg] - path to output bam
    :param mmlim: [char/arg] - mismatch limit
    :param bLen: [int8_t/arg] - barcode length
    :return: [cystr]
    """
    return BamRescue(inBam, outBam, mmlim, bLen)


cdef cystr BamRescue(cystr inBam,
                     cystr outBam,
                     char mmlim, int8_t bLen):
    input_bam = pysam.AlignmentFile(inBam, "rb")
    output_bam = pysam.AlignmentFile(outBam, "wb", template=input_bam)
    cdef AlignedSegment_t read
    cdef list recList
    cdef int RefID, RNext, Pos, MPos
    cdef bint IsRead1, IsRev, Pass
    obw = output_bam.write
    for Pass, gen in groupby(input_bam, SKIP_READS):
        if not Pass:
            [obw(read) for read in gen]
            continue
        for RefID, gen1 in groupby(input_bam, REF_ID):
            for Pos, gen2 in groupby(gen1, POS):
                for RNext, gen3 in groupby(gen2, RNEXT):
                    for MPos, gen4 in groupby(gen3, MPOS):
                        for IsRead1, gen5 in groupby(gen4, IS_READ1):
                            for IsRev, FinalGen in groupby(gen5, IS_REV):
                                recList = list(FinalGen)
                                if len(recList) == 1:
                                    obw(recList[0])
                                    continue
                                recList = BFF(recList, bLen, mmlim)
                                [obw(read) for read in recList]
    input_bam.close()
    output_bam.close()
    return output_bam.filename


cdef inline int8_t StringHD(char * str1, char * str2, int8_t bLen) nogil:
    cdef size_t index
    cdef int8_t ret = 0
    for index in range(bLen):
        if(str1[index] != str2[index]):
            ret += 1
    return ret


cdef inline pBarcodeHD(AlignedSegment_t query, AlignedSegment_t cmp,
                       int8_t bLen):
    cdef cystr BC1, BC2
    BC1 = query.query_name
    BC2 = cmp.query_name
    return StringHD(<char*> BC1, <char*> BC2, bLen)

cdef AlignedSegment_t MergeBamRecList(list recList):
    return recList[0]

'''
cdef inline pBarcodeHD(AlignedSegment_t query, AlignedSegment_t cmp,
                       int8_t bLen):
    return BarcodeHD(query._delegate, cmp._delegate, bLen)
'''
