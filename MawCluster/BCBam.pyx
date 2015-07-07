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
import pysam
import cython

from .BCFastq import GetDescriptionTagDict as getdesc
from . import BCFastq
from utilBMF.HTSUtils import (printlog as pl,
                              FractionAligned, FractionSoftClipped,
                              SWRealignAS, pPileupRead, BedtoolsBam2Fq,
                              BwaswCall, samtoolsMergeBam, pFastqProxy,
                              TrimExt, align_bwa_mem)
from utilBMF.ErrorHandling import (IllegalArgumentError, ThisIsMadness as Tim,
                                   MissingExternalTool)
from .SVUtils import returnDefault
from utilBMF import HTSUtils
from warnings import warn
import SecC

cimport numpy as np
cimport cython

cimport utilBMF.HTSUtils
ctypedef utilBMF.HTSUtils.pFastqProxy pFastqProxy_t
npchararray = char.array
oagseq = oag("seq")
oagqqual = oag("query_qualities")
mcoptBS = mc("opt", "BS")
# BSStringList = map(mcoptBS, readList)


def AbraCadabra(inBAM, outBAM="default",
                jar="default", memStr="default", ref="default",
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
    if(memStr == "default"):
        memStr = "-Xmx16G"
        pl("Default memory string used: " + memStr)
    else:
        pl("Non-default memory string used: " + memStr)
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
        "memStr=\"{}\", ref=\"{}\", threads=\"{}\", ".format(memStr,
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
    command = ("java {} -jar {} --in {}".format(memStr, jar, inBAM) +
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


def BarcodeSort(inBAM, outBAM="default", paired=True):
    if(outBAM == "default"):
        outBAM = '.'.join(inBAM.split('.')[0:-1]) + "barcodeSorted.bam"
    pl("BarcodeSort. Input: {}. Output: {}.".format(inBAM, outBAM))
    outsam = '.'.join(outBAM.split('.')[0:-1]) + ".sam"
    headerCommand = "samtools view -H {}".format(inBAM)
    pl(headerCommand)
    check_call(shlex.split(headerCommand), shell=False, stdout=outsam)
    pl("Now converting bam to sam for sorting by barcode.")
    if(paired is False):
        cmd = ("samtools view {} | ".format(inBAM) +
               "awk 'BEGIN {{FS=\"\t\";OFS=\"\t\"}};{{print "
               "$(NF-2),$0}}' - | sort | cut -f2- >> {}".format(outsam))
    return outBAM


def mergeBarcodes(reads1, reads2, outfile="default"):
    cdef cystr concatBarcode
    cdef AlignedSegment_t entry1, entry2
    pl("mergeBarcodes. R1: {}. R2: {}".format(reads1, reads2))
    reader1 = pysam.AlignmentFile(reads1, "rb")
    reader2 = pysam.AlignmentFile(reads2, "rb")
    if(outfile == "default"):
        outfile = '.'.join(reads1.split('.')[0:-2]) + '.merged.bam'
    outSAM = pysam.AlignmentFile(outfile, "wb", template=reader1)
    osw = outSAM.write
    for entry1 in reader1:
        entry2 = reader2.next()
        assert entry1.qname == entry2.qname
        # print("Barcode 1: {}. Barcode 2: {}.".format(Barcode1,Barcode2))
        concatBarcode = entry1.opt("BS") + entry2.opt("BS")
        # print("New barcode will be {}".format(concatBarcode))
        entry1.setTag("BS", concatBarcode)
        entry2.setTag("BS", concatBarcode)
        osw(entry1)
        osw(entry2)
    reader1.close()
    reader2.close()
    outSAM.close()
    return outfile


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


def pairedBarcodeTagging(
        cystr fq1,
        cystr fq2,
        cystr bam,
        cystr outBAMFile="default",
        cystr suppBam="default",
        cystr conversionXml="default", cystr realigner="none",
        cystr ref="default"):
    """
    TODO: Unit test for this function.
    Tags a BAM from comments in fastq file records.
    Not used for reads where the comment has
    been appended to each record in the sam file (e.g., bwa mem -C)
    :param cystr fq1 - Path to read 1 fastq
    :param cystr fq2 - Path to read 2 fastq
    :param cystr bam - Path to aligned, name-sorted bam
    :param cystr conversionXml - Path to xml from bcl2fastq for
    annotating reads.
    :param cystr realigner - if "gatk" in this string, add RG.
    """
    cdef ndarray[np.int64_t, ndim=1] PhredQuals1, PhredQuals2, FA1, FA2
    cdef AlignedSegment_t entry, read1bam, read2bam
    cdef double r1FracAlign, r2FracAlign, r1FracSC, r2FracSC
    cdef int FM, ND1, ND2
    cdef bint addDefault, bwaswRescue, passing
    cdef cystr coorString, cStr, contigSetStr
    cdef dict descDict1, descDict2
    cdef pFastqProxy_t pFq1, pFq2
    # cdef AlignmentFile postFilterBAM, outBAM, suppBAM
    if(outBAMFile == "default"):
        outBAMFile = '.'.join(bam.split('.')[0:-1]) + ".tagged.bam"
    if(suppBam == "default"):
        suppBam = '.'.join(bam.split('.')[0:-1]) + '.2ndSupp.bam'
    pl("pairedBarcodeTagging. Input bam: %s. outputBAM: %s" % (bam,
                                                               outBAMFile))
    cStr = "pairedBarcodeTagging({}, {}, {})".format(fq1, fq2,
                                                     bam)
    pl("Command string to reproduce call: {}".format(cStr))
    pl("realigner: %s" % realigner)
    # read1Handle = pysam.FastqFile(fq1)
    # read2Handle = pysam.FastqFile(fq2)
    read1Handle = pysam.FastqFile(fq1)
    read2Handle = pysam.FastqFile(fq2)
    postFilterBAM = pysam.AlignmentFile(bam, "rb")
    outBAM = pysam.AlignmentFile(outBAMFile, "wb", template=postFilterBAM)
    suppBAM = pysam.AlignmentFile(suppBam, "wb", template=postFilterBAM)
    if(conversionXml != "default"):
        convData = SecC.SecC.BuildRunDict(conversionXml)
    obw = outBAM.write
    addDefault = ("gatk" in realigner)
    r1hn = read1Handle.next
    r2hn = read2Handle.next
    for entry in postFilterBAM:
        if(entry.is_secondary or entry.flag >= 2048):
            suppBAM.write(entry)
            continue
        if(entry.is_read1):
            read1bam = entry
            pFq1 = pFastqProxy.fromFastqProxy(r1hn())
            continue
        elif(entry.is_read2):
            read2bam = entry
            pFq2 = pFastqProxy.fromFastqProxy(r2hn())
        descDict1 = getdesc(pFq1.comment)
        descDict2 = getdesc(pFq2.comment)
        passing = "Pass" in descDict1["FP"]
        FM = int(descDict1["FM"])
        try:
            ND1 = int(descDict1["ND"])
            ND2 = int(descDict2["ND"])
            PhredQuals1 = np.array(descDict1["PV"].split(","), dtype=np.int64)
            PhredQuals2 = np.array(descDict2["PV"].split(","), dtype=np.int64)
            FA1 = np.array(descDict1["FA"].split(","), dtype=np.int64)
            FA2 = np.array(descDict2["FA"].split(","), dtype=np.int64)
        except KeyError:
            raise Tim("Number of Differences tag required for "
                      "BMFTools >= v0.0.7")
        except ValueError:
            raise ValueError("ND tag value is invalid: "
                             "%s %s" % (descDict1["ND"], descDict2["ND"]))
        # If the read is reversed, the PV tag must be reversed to match
        if(read1bam.is_reverse):
            PhredQuals1 = PhredQuals1[::-1]
            FA1 = FA1[::-1]
        if(read2bam.is_reverse):
            PhredQuals2 = PhredQuals2[::-1]
            FA2 = FA2[::-1]
        r1FracAlign = FractionAligned(read1bam)
        r1FracSC = FractionSoftClipped(read1bam)
        """
        if(r1FracAlign < minAF and not read1bam.is_unmapped):
            read1bam = SWRealignAS(read1bam, postFilterBAM, ref=ref)
            r1FracAlign = FractionAligned(read1bam)
        """
        r2FracAlign = FractionAligned(read2bam)
        r2FracSC = FractionSoftClipped(read2bam)
        """
        if(r2FracAlign < minAF and not read2bam.is_unmapped):
            read2bam = SWRealignAS(read2bam, postFilterBAM, ref=ref)
            r2FracAlign = FractionAligned(read2bam)
        """
        coorString = ",".join(sorted([":".join([PysamToChrDict[
            read1bam.reference_id], str(read1bam.pos)]), ":".join([
                PysamToChrDict[read2bam.reference_id], str(read2bam.pos)])]))
        contigSetStr = ",".join(sorted(
            [PysamToChrDict[read1bam.reference_id],
             PysamToChrDict[read2bam.reference_id]]))

        if(addDefault):
            read1bam.set_tags([("RP", coorString, "Z"),
                               ("SC", contigSetStr, "Z"),
                               ("FM", FM, "i"),
                               ("BS", descDict1["BS"], "Z"),
                               ("FP", passing, "i"),
                               ("PV", ",".join(PhredQuals1.astype(str)), "Z"),
                               ("FA", ",".join(FA1.astype(str)), "Z"),
                               ("ND", ND1, "i"),
                               ("NF", ND1 * 1. / FM, "f"),
                               ("RG", "default", "Z"),
                               ("AF", r1FracAlign, "f"),
                               ("SF", r1FracSC, "f")])
            read2bam.set_tags([("RP", coorString, "Z"),
                               ("SC", contigSetStr, "Z"),
                               ("FM", FM, "i"),
                               ("BS", descDict1["BS"], "Z"),
                               ("FP", passing, "i"),
                               ("PV", ",".join(PhredQuals2.astype(str)), "Z"),
                               ("FA", ",".join(FA2.astype(str)), "Z"),
                               ("ND", ND2, "i"),
                               ("NF", ND2 * 1. / float(FM), "f"),
                               ("RG", "default", "Z"),
                               ("AF", r2FracAlign, "f"),
                               ("SF", r2FracSC, "f")])
        else:
            read1bam.set_tags([("RP", coorString, "Z"),
                               ("SC", contigSetStr, "Z"),
                               ("FM", FM, "i"),
                               ("BS", descDict1["BS"], "Z"),
                               ("FP", int("Pass" in descDict1["FP"]), "i"),
                               ("PV", ",".join(PhredQuals1.astype(str)), "Z"),
                               ("FA", ",".join(FA1.astype(str)), "Z"),
                               ("ND", ND1, "i"),
                               ("NF", 1. * ND1 / FM, "f"),
                               ("AF", r1FracAlign, "f"),
                               ("SF", r1FracSC, "f")])
            read2bam.set_tags([("RP", coorString, "Z"),
                               ("SC", contigSetStr, "Z"),
                               ("FM", FM, "i"),
                               ("BS", descDict1["BS"], "Z"),
                               ("FP", int("Pass" in descDict1["FP"]), "i"),
                               ("PV", ",".join(PhredQuals2.astype(str)), "Z"),
                               ("FA", ",".join(FA2.astype(str)), "Z"),
                               ("ND", ND2, "i"),
                               ("NF", 1. * ND2 / FM, "f"),
                               ("AF", r2FracAlign, "f"),
                               ("SF", r2FracSC, "f")])
        read1bam.is_qcfail = (not passing)
        read2bam.is_qcfail = (not passing)
        # I used to mark the BAMs at this stage, but it's not appropriate to
        # do so until after indel realignment.
        obw(read1bam)
        obw(read2bam)
    suppBAM.close()
    outBAM.close()
    postFilterBAM.close()
    return outBAMFile


def compareRecs(RecordList, oagseq=oagseq, oagqqual=oagqqual):
    Success = True
    seqs = map(oagseq, RecordList)
    seqs = [str(record.seq) for record in RecordList]
    stackArrays = tuple([npchararray(s, itemsize=1) for s in seqs])
    seqArray = nvstack(stackArrays)
    # print(repr(seqArray))

    quals = np.array(map(oagqqual, RecordList))
    qualA = ccopy(quals)
    qualC = ccopy(quals)
    qualG = ccopy(quals)
    qualT = ccopy(quals)
    qualA[seqArray != "A"] = 0
    qualASum = nsum(qualA, 0)
    qualC[seqArray != "C"] = 0
    qualCSum = nsum(qualC, 0)
    qualG[seqArray != "G"] = 0
    qualGSum = nsum(qualG, 0)
    qualT[seqArray != "T"] = 0
    qualTSum = nsum(qualT, 0)
    qualAllSum = nvstack([qualASum, qualCSum, qualGSum, qualTSum])
    newSeq = "".join(map(Num2Nuc, nargmax(qualAllSum, 0)))
    MaxPhredSum = np.amax(qualAllSum, 0)  # Avoid calculating twice.
    phredQuals = nsub(nmul(2, MaxPhredSum),
                      nsum(qualAllSum, 0))
    phredQuals[phredQuals < 0] = 0
    outRec = RecordList[0]
    outRec.seq = newSeq
    if(np.any(np.greater(phredQuals, 93))):
        outRec.setTag("PV", ",".join(phredQuals.astype(str)))
    phredQuals[phredQuals > 93] = 93
    outRec.query_qualities = phredQuals
    return outRec, Success


def ConsolidateInferred(inBAM, outBAM="default"):
    if(outBAM == "default"):
        outBAM = '.'.join(inBAM.split('.')[0:-1]) + "consolidated.bam"
    inputHandle = pysam.AlignmentFile(inBAM, 'rb')
    outputHandle = pysam.AlignmentFile(outBAM, 'wb', template=inputHandle)
    workBC1 = ""
    workBC2 = ""
    Set1 = []
    Set2 = []
    for record in inputHandle:
        if(record.is_read1):
            barcodeRecord1 = record.opt("RP")
            if(workBC1 == ""):
                workBC1 = barcodeRecord1
                Set1 = []
                Set1.append(record)
            elif(workBC1 == barcodeRecord1):
                Set1.append(record)
            else:
                mergeRec1, success = compareRecs(Set1)
                if(success is False):
                    mergeRec1.setTag("FP", 0)
                outputHandle.write(mergeRec1)
                Set1 = [record]
                workBC1 = barcodeRecord1
        if(record.is_read2):
            barcodeRecord2 = record.opt("RP")
            if(workBC2 == ""):
                workBC2 = barcodeRecord2
                Set2 = []
                Set2.append(record)
            elif(workBC2 == barcodeRecord2):
                Set2.append(record)
            else:
                mergeRec2, success = compareRecs(Set2)
                if(success is False):
                    mergeRec2.setTag("FP", 0)
                outputHandle.write(mergeRec2)
                Set2 = [record]
                workBC2 = barcodeRecord2
    inputHandle.close()
    outputHandle.close()
    return outBAM


def singleBarcodeTagging(cystr fastq, cystr bam, cystr outputBAM="default",
                         cystr suppBam="default"):
    cdef pFastqProxy_t FqPrx
    cdef pysam.cfaidx.FastqProxy tempRead
    cdef AlignedSegment_t entry
    cdef pysam.cfaidx.FastqFile reads
    cdef dict descDict
    """
    TODO: Unit test for this function.
    """
    if(outputBAM == "default"):
        outputBAM = TrimExt(bam) + ".tagged.bam"
    if(suppBam == "default"):
        suppBam = TrimExt(bam) + '.2ndSupp.bam'
    pl("singleBarcodeTagging. Fq: {}. outputBAM: {}".format(bam, outputBAM))
    reads = pysam.FastqFile(fastq)
    postFilterBAM = pysam.AlignmentFile(bam, "rb")
    suppBAM = pysam.AlignmentFile(suppBam, "wb", template=postFilterBAM)
    outBAM = pysam.AlignmentFile(outputBAM, "wb", template=postFilterBAM)
    for entry in postFilterBAM:
        if(entry.is_secondary or entry.flag & 2048):
            suppBAM.write(entry)
            continue
        else:
            try:
                tempRead = reads.next()
                FqPrx = pFastqProxy.fromFastqProxy(tempRead)
            except StopIteration:
                break
        descDict = getdesc(FqPrx.comment)
        for key in descDict.iterkeys():
            entry.setTag(key, descDict[key])
        if("Pass" in descDict["FP"]):
            entry.set_tag("FP", 1)
        else:
            entry.set_tag("FP", 0)
            entry.is_qcfail = True
        outBAM.write(entry)
    outBAM.close()
    postFilterBAM.close()
    return outputBAM


@cython.returns(bint)
def FracSoftclippedTest(AlignedSegment_t rec,
                        double maxFracSoftClipped=0.25):
    if(FractionSoftClipped(rec) >= maxFracSoftClipped):
        return False
    return True


def GetFracSCPartial(double maxFracSoftClipped):
    """
    Returns a partial for FracSoftclippedTest so that it can
    be passed into AbstractBamFilter.
    """
    return partial(FracSoftclippedTest,
                   maxFracSoftClipped=maxFracSoftClipped)


def AbstractBamFilter(inBAM, failBAM="default", passBAM="default",
                      func=returnDefault, appendStr=""):
    cdef AlignedSegment_t rec, r1, r2
    cdef AlignmentFile inHandle, raHandle, nrHandle
    if(failBAM == "default"):
        failBAM = ".".join(inBAM.split(".")[:-1] + [appendStr, "Fail", "bam"])
    if(passBAM == "default"):
        passBAM = ".".join(inBAM.split(".")[:-1] + [appendStr, "Pass", "bam"])
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    raHandle = pysam.AlignmentFile(failBAM, "wb", template=inHandle)
    nrHandle = pysam.AlignmentFile(passBAM, "wb", template=inHandle)
    pl("Got all handles. Now filtering!")
    for rec in inHandle:
        if(rec.is_read1):
            r1 = rec
            continue
        elif(rec.is_read2):
            r2 = rec
        if(func(r1) or func(r2)):
            raHandle.write(r1)
            raHandle.write(r2)
        else:
            nrHandle.write(r1)
            nrHandle.write(r2)
    return passBAM, failBAM


def GetSoftClips(inBAM, failBAM="default", passBAM="default",
                 double maxFracSoftClipped=0.5):
    """
    Uses the AbstractBamFilter to get reads with Softclipped Fraction >= 0.25
    """
    return AbstractBamFilter(inBAM, failBAM=failBAM, passBAM=passBAM,
                             func=GetFracSCPartial(maxFracSoftClipped),
                             appendStr="SFlt%s" % maxFracSoftClipped)


def AddRATag(inBAM, inplace=False, outBAM="default", RATag="bwasw"):
    """
    Uses sed and samtools view to append a tag to each line of the file
    not in the header.
    """
    tmpfile = str(uuid.uuid4().get_hex()[0:8]) + '.bam'
    tag = "RA:Z:" + RATag
    if(inplace):
        pl("Adding RA Tag 'in-place'.")
    else:
        if(outBAM == "default"):
            outBAM = ".".join(outBAM.split("."))[:-1] + "."
    pl("Adding RA:z:bwasw tag.")
    cStr = ("samtools view -h %s | " % inBAM +
            "awk 'FS=OFS=\"\t\" {{if($1 !~ \"^@\") {{print $0, "
            "\"RA:Z:bwasw\"}} else {{print $0}}}}'"
            " | samtools view -Sbh - > %s" % (tag, tmpfile))

    check_call(cStr, shell=True)
    if(inplace):
        check_call(["mv", tmpfile, inBAM])
        return inBAM
    else:
        check_call(["mv", tmpfile, outBAM])
        return outBAM


def RealignSFReads(inBAM, double maxFracSoftClipped=0.5,
                   ref="default", outBAM="default"):
    """
    Realigns reads which have a Softclipped Fraction that is above
    maxFracSoftClipped
    """
    if(outBAM == "default"):
        outBAM = ".".join(inBAM.split()[:-1]) + ".SWRealigned.bam"
    if(ref == "default"):
        raise Tim("ref required for bwasw alignment.")
    print("Getting soft-clipped reads for bwasw realignment")
    NoRealign, Realign = GetSoftClips(
        inBAM, maxFracSoftClipped=maxFracSoftClipped)
    pl("Now converting bam to fastq")
    ReadFastq1, ReadFastq2 = BedtoolsBam2Fq(Realign)
    pl("bwasw call!")
    RealignedBAM = BwaswCall(ReadFastq1, ReadFastq2, ref=ref)
    pl("Sorting the unrealigned bam!")
    SortNoRealign = HTSUtils.CoorSortAndIndexBam(NoRealign, delete=True)
    pl("Sorting the realigned bam!")
    SortRealign = HTSUtils.CoorSortAndIndexBam(Realign, delete=True)
    pl("Merging the bams!")
    samtoolsMergeBam([SortNoRealign, SortRealign],
                     outBAM=outBAM)
    pl("Adding the RA tag")
    AddRATag(outBAM, inplace=True, RATag="bwasw")
    return outBAM


cdef dict cGetCOTagDict(AlignedSegment_t read):
    cdef cystr s, cStr
    cStr = read.opt("CO")
    return dict([s.split("=") for s in cStr.split("|")[1:]])


cpdef dict pGetCOTagDict(AlignedSegment_t read):
    return cGetCOTagDict(read)


cpdef AlignedSegment_t pTagAlignedSegmentHG37(AlignedSegment_t read):
    return TagAlignedSegmentHG37(read)

'''
@cython.boundscheck(False)
@cython.wraparound(False)
'''
cdef cystr RPStringNonHG37(AlignedSegment_t read, dict RefIDDict=None):
    return (RefIDDict[read.reference_id] + ":%s," % read.pos +
            RefIDDict[read.rnext] +
            ":%s" % read.mpos)


'''
@cython.boundscheck(False)
@cython.wraparound(False)
'''
cdef AlignedSegment_t TagAlignedSegment(
        AlignedSegment_t read, dict RefIDDict=None):
    """
    Adds necessary information from a CO: tag to appropriate other tags.
    """
    assert RefIDDict is not None
    cdef dict CommentDict
    cdef int FM, ND, FPInt
    cdef double NF, AF, SF
    cdef ndarray[np.int64_t, ndim=1] PhredQuals, FA
    CommentDict = cGetCOTagDict(read)
    PhredQuals = np.array(CommentDict["PV"].split(","), dtype=np.int64)
    FA = np.array(CommentDict["FA"].split(","), dtype=np.int64)
    if(read.is_reverse):
        PhredQuals = PhredQuals[::-1]
        FA = FA[::-1]
    ND = int(CommentDict["ND"])
    FM = int(CommentDict["FM"])
    FPInt = 1 if("Pass" in CommentDict["FP"]) else 0
    if(not FPInt):
        read.flag += 512
    NF = ND * 1. / FM
    AF = getAF(read)
    SF = getSF(read)
    read.set_tags([("BS", CommentDict["BS"], "Z"),
                   ("FM", FM, "i"),
                   ("PV", ",".join(PhredQuals.astype(str)), "Z"),
                   ("FA", ",".join(FA.astype(str)), "Z"),
                   ("FP", FPInt, "i"),
                   ("ND", int(CommentDict["ND"]), "i"),
                   ("NF", NF, "f"),
                   ("AF", AF, "f"),
                   ("SF", SF, "f"),
                   ("RP", RPStringNonHG37(read, RefIDDict), "Z")
                   ] + read.get_tags())
    read.set_tag("CO", None)  # Delete the CO tag.
    return read


cdef AlignedSegment_t TagAlignedSegmentHG37(
        AlignedSegment_t read):
    """
    Adds necessary information from a CO: tag to appropriate other tags.
    """
    cdef dict CommentDict
    cdef int FM, ND, FPInt
    cdef double NF, AF, SF
    cdef ndarray[np.int64_t, ndim=1] PhredQuals, FA
    CommentDict = cGetCOTagDict(read)
    PhredQuals = np.array(CommentDict["PV"].split(","), dtype=np.int64)
    FA = np.array(CommentDict["FA"].split(","), dtype=np.int64)
    if(read.is_reverse):
        PhredQuals = PhredQuals[::-1]
        FA = FA[::-1]
    ND = int(CommentDict["ND"])
    FM = int(CommentDict["FM"])
    FPInt = 1 if("Pass" in CommentDict["FP"]) else 0
    if(not FPInt):
        read.flag += 512
    NF = ND * 1. / FM
    AF = getAF(read)
    SF = getSF(read)
    read.set_tags([("BS", CommentDict["BS"], "Z"),
                   ("FM", FM, "i"),
                   ("PV", ",".join(PhredQuals.astype(str)), "Z"),
                   ("FA", ",".join(FA.astype(str)), "Z"),
                   ("FP", FPInt, "i"),
                   ("ND", int(CommentDict["ND"]), "i"),
                   ("NF", NF, "f"),
                   ("AF", AF, "f"),
                   ("SF", SF, "f"),
                   ("RP", RPString(read), "Z")
                   ] + read.get_tags())
    read.set_tag("CO", None)  # Delete the CO tag.
    return read


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


cdef class TagBamPipe:
    """
    Doesn't require a precomputed dictionary, since it builds
    it from the input stream.
    """
    def __init__(self, bint bin_input, bint bin_output,
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
        self.RefIDDict = dict(list(enumerate(self.inHandle.references)) +
                              [(-1, "*")])

    cdef write(self, AlignedSegment_t read):
        self.outHandle.write(read)

    cpdef process(self):
        cdef AlignedSegment_t read
        for read in self.inHandle:
            read = TagAlignedSegment(read, RefIDDict=self.RefIDDict)
            self.write(read)

cdef class TagBamPipeHG37(BamPipe):
    """
    Not unlike TagBamPipe, but memoizes HG37's ref id dictionary and
    inlines it into switches.
    """
    def __init__(self, bint bin_input, bint bin_output,
                 bint uncompressed_output=False):
        super(BamPipe, self).__init__(TagAlignedSegmentHG37,
                                      bin_input, bin_output,
                                      uncompressed_output=uncompressed_output)


def PipeBarcodeTagCOBam(char flag):
    """
    Takes a SAM input stream, tags the bam, and converts
    it into a bam, all in one fell swoop!
    Uses a bitwise flag.
    flag & 4 is true to emit uncompressed
    flag & 2 is true if input is textual (sam) format
    flag & 1 is true if output is textual (sam) format
    """
    from sys import stderr
    cdef AlignmentFile_t inHandle, outHandle
    cdef AlignedSegment_t read
    cdef dict RefIDDict
    cdef cystr input_mode
    input_mode = "r" if(flag & 2) else "rb"
    if(flag & 1):
        output_mode = "w"
    elif(flag & 4):
        output_mode = "wbu"
    else:
        output_mode = "wb"
    stderr.write("Now tagging a piped BAM with input mode "
                 "%s and output mode %s\n" % (input_mode, output_mode))
    inHandle = pysam.AlignmentFile("-", input_mode)
    outHandle = pysam.AlignmentFile("-", output_mode, template=inHandle)
    RefIDDict = dict(list(enumerate(inHandle.references)) + [(-1, "*")])
    [outHandle.write(TagAlignedSegment(read, RefIDDict)) for read in inHandle]
    inHandle.close()
    outHandle.close()
    return 1


cdef cystr cBarcodeTagCOBam(AlignmentFile inbam,
                            AlignmentFile outbam):
    """In progress
    """
    cdef AlignedSegment_t read
    cdef dict RefIDDict = dict(list(enumerate(inbam.references)) +
                               [(-1, "*")])
    for read in inbam:
        outbam.write(TagAlignedSegment(read, RefIDDict))
    inbam.close()
    outbam.close()
    return outbam.filename


cpdef cystr pBarcodeTagCOBam(cystr bam, cystr outbam=None):
    """In progress
    """
    pl("Tagging BAM with barcode tags in CO field ...")
    cdef AlignmentFile inHandle
    inHandle = pysam.AlignmentFile(bam, "rb")
    outbam = ".".join(bam.split(".")[:-1]) + ".tagged.bam"
    return cBarcodeTagCOBam(
        inHandle, pysam.AlignmentFile(outbam, "wb", template=inHandle))


def AlignAndTagMem(*args, **kwargs):
    raise DeprecationWarning("AlignAndTagMem is deprecated. "
                             "Please use utilBMF.HTSUtils.PipeAlignTag.")


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


cdef cystr BamRescue(cystr inBam,
                     cystr outBam,
                     char mmlim, int8_t bLen):
    input_bam = pysam.AlignmentFile(inBam, "rb")
    output_bam = pysam.AlignmentFile(outBam, "wb", template=input_bam)
    cdef AlignedSegment_t query_read, cmp_read, read
    cdef ndarray[char, ndim=2] distmtx
    cdef py_array arr, nr_arr
    cdef list recList, merging_arr_sets
    cdef size_t size, qctr, cmpctr, t
    cdef object gen
    cdef int8_t dist
    cdef set rsqidx
    cdef int RefID, RNext, Pos, MPos
    cdef bint IsRead1
    cdef object cfi = chain.from_iterable
    obw = output_bam.write
    for RefID, gen in groupby(input_bam, REF_ID):
        for Pos, gen1 in groupby(gen, POS):
            for RNext, gen2 in groupby(gen1, RNEXT):
                for MPos, gen3 in groupby(gen2, MPOS):
                    for IsRead1, FinalGen in groupby(gen3, IS_READ1):
                        setsToMerge = {}
                        recList = list(gen)
                        size = len(recList)
                        distmtx = np.zeros([size, size], dtype=np.int8)
                        cmpctr = 1
                        for qctr, query_read in enumerate(recList):
                            for cmp_read in recList[qctr + 1:]:
                                distmtx[qctr][cmpctr] = pBarcodeHD(
                                        query_read, cmp_read, bLen)
                                cmpctr += 1
                        cmpctr = 1
                        merging_arr_sets = []
                        rsqidx = set([])
                        ria = rsqidx.add
                        for qctr in range(size):
                            for cmpctr in range(qctr + 1):
                                if(distmtx[qctr][cmpctr] < mmlim):
                                    arr = array('i', [
                                        t for t, dist in
                                        enumerate(distmtx[qctr][:]) if
                                        dist < mmlim and
                                        t not in rsqidx])
                                    merging_arr_sets.append(arr)
                                    [ria(t) for t in arr]
                        rsqidx = array('i', list(cfi(merging_arr_sets)))
                        nr_arr = array('i', [t for t in range(size) if
                                             t not in rsqidx])
                        [obw(recList[t]) for t in range(size) if
                         t not in rsqidx]
                        for arr in merging_arr_sets:
                            read = BarcodeRescueBam(recList[t] for t in arr)
                            obw(read)
    return output_bam.filename


cdef BarcodeRescueBam(list recList):
    return recList[0]


cdef pBarcodeHD(AlignedSegment_t query, AlignedSegment_t cmp, int8_t bLen):
    cdef bam1_t * query_src
    cdef bam1_t * cmp_src
    query_src = query._delegate
    cmp_src = cmp._delegate
    return BarcodeHD(query_src, cmp_src, bLen)
