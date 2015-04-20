#!python
# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True
from __future__ import division
import shlex
import subprocess
import os
import shutil
import logging
from copy import copy as ccopy
from os import path
import operator
from operator import attrgetter as oag
from operator import div as odiv
import string
import uuid
import sys

import numpy as np
from numpy import array as nparray
from numpy import sum as nsum
from numpy import multiply as nmultiply
from numpy import subtract as nsubtract
from numpy import argmax as nargmax
from numpy import vstack as nvstack
from numpy import char
from Bio import SeqIO
import pysam
import cython
from cytoolz import map as cmap

from MawCluster.BCFastq import letterNumDict, GetDescriptionTagDict as getdesc
from MawCluster import BCFastq
from utilBMF.HTSUtils import (printlog as pl, PysamToChrDict, ThisIsMadness,
                              FractionAligned, FractionSoftClipped,
                              SWRealignAS, pPileupRead)
from utilBMF.ErrorHandling import IllegalArgumentError
from utilBMF import HTSUtils
import SecC
cimport numpy as np
cimport cython
cimport pysam.calignmentfile
ctypedef pysam.calignmentfile.AlignedSegment cAlignedSegment
npchararray = char.array


@cython.locals(fixMate=cython.bint)
def AbraCadabra(inBAM, outBAM="default",
                jar="default", memStr="default", ref="default",
                threads="4", bed="default", working="default",
                log="default", fixMate=True, tempPrefix="tmpPref",
                rLen=-1, intelPath="default"):
    """
    Calls abra for indel realignment. It supposedly
    out-performs GATK's IndelRealigner.
    Note: bed file must be first 3 columns only and
    coordinate sorted. You will likely need an additional bed file for this.
    """
    if(rLen < 0):
        raise IllegalArgumentError("Read length must be set to call abra due"
                                   " to the benefits of inferring ideal para"
                                   "meters from the !")
    if(jar == "default"):
        HTSUtils.FacePalm("Required: Path to abra jar!")
    else:
        pl("Non-default abra jar used: " + jar)
    if(memStr == "default"):
        memStr = "-Xmx16G"
        pl("Default memory string used: " + memStr)
    else:
        pl("Non-default memory string used: " + memStr)
    if(ref == "default"):
        raise ValueError("Reference fasta must be provided!")
    else:
        pl("Reference file set: {}.".format(ref))
    if(bed == "default"):
        raise ValueError("Bed file required.")
    else:
        pl("Bed file set: {}.".format(bed))
    if(working == "default"):
        working = inBAM.split('.')[0] + ".working_dir"
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
    bed = AbraKmerBedfile(bed, ref=ref, abra=jar,
                          rLen=rLen)
    if(path.isfile(inBAM + ".bai") is False):
        pl("No bam index found for input bam - attempting to create.")
        subprocess.check_call(['samtools', 'index', inBAM])
        if(path.isfile(inBAM + ".bai") is False):
            inBAM = HTSUtils.CoorSortAndIndexBam(inBAM, outBAM, uuid=True)
    command = ("java {} -jar {} --in {}".format(memStr, jar, inBAM) +
               " --out {} --ref {} --targets".format(outBAM, ref) +
               " {} --threads {} ".format(bed, threads) +
               "--working {} --mbq 200".format(working))
    pl("Command: {}.".format(command))
    subprocess.check_call(shlex.split(command), shell=False)
    pl("Deleting abra's intermediate directory.")
    subprocess.check_call(["rm", "-rf", working])
    if(fixMate):
        pl("Now fixing mates after abra's realignment.")
        tempFilename = tempPrefix + str(
            uuid.uuid4().get_hex()[0:8]) + ".working.tmp"
        nameSorted = HTSUtils.NameSort(outBAM)
        commandStrFM = "samtools fixmate %s %s -O bam" % (nameSorted,
                                                          tempFilename)
        subprocess.check_call(shlex.split(commandStrFM))
        subprocess.check_call(["rm", "-rf", nameSorted])
        subprocess.check_call(["mv", tempFilename, outBAM])
    return outBAM


@cython.locals(rLen=cython.long)
def AbraKmerBedfile(inbed, rLen=-1, ref="default", outbed="default",
                    nt=4, abra="default"):
    if(abra == "default"):
        raise ThisIsMadness(
            "Path to abra jar required for running KmerSizeCalculator.")
    if(ref == "default"):
        raise ThisIsMadness(
            "Path to reference required for running KmerSizeCalculator.")
    if(inbed == "default"):
        raise ThisIsMadness(
            "Path to input bed required for running KmerSizeCalculator.")
    if(rLen < 0):
        raise ThisIsMadness(
            "Read length required for running KmerSizeCalculator.")
    if(outbed == "default"):
        outbed = ".".join(inbed.split(".")[0:-1] + ["abra", "kmer"])
    commandStr = ("java -cp %s abra.KmerSizeEvaluator " % abra +
                  "%s %s %s %s %s" % (rLen, ref, outbed, nt, inbed))
    pl("AbraKmerSizeEvaluator call string: %s" % commandStr)
    subprocess.check_call(shlex.split(commandStr))
    return outbed


def Bam2Sam(inBAM, outsam):
    pl("Bam2Sam. Input: {}. Output: {}.".format(inBAM, outsam))
    output = open(outsam, 'w', 0)
    command_str = 'samtools view -h {}'.format(inBAM)
    pl(command_str)
    subprocess.check_call(shlex.split(command_str), stdout=output, shell=False)
    return(command_str, outsam)


def BarcodeSort(inBAM, outBAM="default", paired=True):
    if(outBAM == "default"):
        outBAM = '.'.join(inBAM.split('.')[0:-1]) + "barcodeSorted.bam"
    pl("BarcodeSort. Input: {}. Output: {}.".format(inBAM, outBAM))
    outsam = '.'.join(outBAM.split('.')[0:-1]) + ".sam"
    headerCommand = "samtools view -H {}".format(inBAM)
    pl(headerCommand)
    subprocess.check_call(shlex.split(headerCommand),
                          shell=False, stdout=outsam)
    pl("Now converting bam to sam for sorting by barcode.")
    if(paired is False):
        cmd = ("samtools view {} | ".format(inBAM) +
               "awk 'BEGIN {{FS=\"\t\";OFS=\"\t\"}};{{print "
               "$(NF-2),$0}}' - | sort | cut -f2- >> {}".format(outsam))
    return outBAM


def mergeBarcodes(reads1, reads2, outfile="default"):
    pl("mergeBarcodes. R1: {}. R2: {}".format(reads1, reads2))
    reader1 = pysam.Samfile(reads1, "rb")
    reader2 = pysam.Samfile(reads2, "rb")
    if(outfile == "default"):
        outfile = '.'.join(reads1.split('.')[0:-2]) + '.merged.bam'
    outSAM = pysam.Samfile(outfile, "wb", template=reader1)
    for entry1 in reader1:
        entry2 = reader2.next()
        assert entry1.qname == entry2.qname
        Barcode1 = entry1.opt("BS")
        Barcode2 = entry1.opt("BS")
        # print("Barcode 1: {}. Barcode 2: {}.".format(Barcode1,Barcode2))
        concatBarcode = Barcode1 + Barcode2
        # print("New barcode will be {}".format(concatBarcode))
        entry1.setTag("BS", concatBarcode)
        entry2.setTag("BS", concatBarcode)
        outSAM.write(entry1)
        outSAM.write(entry2)
    reader1.close()
    reader2.close()
    outSAM.close()
    return outfile


def GATKIndelRealignment(inBAM, gatk="default", ref="default",
                         bed="default", dbsnp="default"):
    if(ref == "default"):
        raise ThisIsMadness("Reference file required for Indel Realignment")
    if(bed == "default"):
        raise ThisIsMadness("Bed file required for Indel Realignment")
    if(gatk == "default"):
        raise ThisIsMadness("Path to GATK Jar required for Indel Realignment")
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
        subprocess.check_call(shlex.split(RTCString))
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
        subprocess.check_call(shlex.split(IRCString))
    except subprocess.CalledProcessError:
        pl("GATK IndelRealignment failed. Still finish the analysis pipeline.")
        return inBAM
    pl("Successful GATK indel realignment. Output: %s" % outBAM)
    return outBAM


def pairedBarcodeTagging(
        cython.str fq1,
        cython.str fq2,
        cython.str bam,
        cython.str outBAMFile="default",
        cython.str suppBam="default",
        cython.str bedfile="default",
        cython.str conversionXml="default", cython.str realigner="default",
        cython.float minAF=0.4, cython.str ref="default"):
    cdef pysam.calignmentfile.AlignedSegment entry, read1bam, read2bam
    cdef cython.float r1FracAlign, r2FracAlign, r1FracSC, r2FracSC
    cdef cython.long FM, ND1, ND2
    cdef cython.bint addDefault, bwaswRescue
    cdef cython.str coorString, cStr, contigSetStr
    cdef dict descDict1, descDict2
    # cdef pysam.calignmentfile.AlignmentFile postFilterBAM, outBAM, suppBAM
    if(realigner == "default"):
        raise ThisIsMadness("realigner must be set to gatk, abra, or none.")
    if(outBAMFile == "default"):
        outBAMFile = '.'.join(bam.split('.')[0:-1]) + ".tagged.bam"
    if(suppBam == "default"):
        suppBam = bam.split('.')[0] + '.2ndSupp.bam'
    pl("pairedBarcodeTagging. Fq: {}. outputBAM: {}".format(bam, outBAMFile))
    cStr = "pairedBarcodeTagging({}, {}, {}, minAF={})".format(fq1, fq2,
                                                               bam, minAF)
    pl("Command string to reproduce call: {}".format(cStr))
    # read1Handle = pysam.FastqFile(fq1)
    # read2Handle = pysam.FastqFile(fq2)
    read1Handle = SeqIO.parse(fq1, "fastq")
    read2Handle = SeqIO.parse(fq2, "fastq")
    # read2Handle = SeqIO.parse(fq2, "fastq")
    postFilterBAM = pysam.AlignmentFile(bam, "rb")
    outBAM = pysam.AlignmentFile(outBAMFile, "wb", template=postFilterBAM)
    suppBAM = pysam.AlignmentFile(suppBam, "wb", template=postFilterBAM)
    if(conversionXml != "default"):
        convData = SecC.SecC.BuildRunDict(conversionXml)
    obw = outBAM.write
    addDefault = (realigner == "gatk")
    r1hn = read1Handle.next
    r2hn = read2Handle.next
    for entry in postFilterBAM:
        if(entry.is_secondary or entry.flag >= 2048):
            suppBAM.write(entry)
            continue
        if(entry.is_read1):
            read1bam = entry
            read1fq = r1hn()
            continue
            # print("Read desc: {}".format(tempRead.description))
        elif(entry.is_read2):
            read2bam = entry
            read2fq = r2hn()
        descDict1 = getdesc(read1fq.description)
        descDict2 = getdesc(read2fq.description)
        FM = int(descDict1["FM"])
        try:
            ND1 = int(descDict1["ND"])
            ND2 = int(descDict2["ND"])
        except KeyError:
            raise ThisIsMadness("Number of Differences tag required for "
                                "BMFTools >= v0.0.7")
        except ValueError:
            raise ValueError("ND tag value is invalid: "
                             "%s %s" % (descDict1["ND"], descDict2["ND"]))
        r1FracAlign = FractionAligned(read1bam)
        r1FracSC = FractionSoftClipped(read1bam)
        if(r1FracAlign < minAF and not read1bam.is_unmapped):
            read1bam = SWRealignAS(read1bam, postFilterBAM, ref=ref)
            r1FracAlign = FractionAligned(read1bam)
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
                               ("FP", int("Pass" in descDict1["FP"]), "i"),
                               ("PV", descDict1["PV"], "Z"),
                               ("FA", descDict1["FA"], "Z"),
                               ("ND", ND1, "i"),
                               ("NF", odiv(ND1, float(FM)), "f"),
                               ("RG", "default", "Z"),
                               ("AF", r1FracAlign, "f"),
                               ("SF", r1FracSC, "f")])
            read2bam.set_tags([("RP", coorString, "Z"),
                               ("SC", contigSetStr, "Z"),
                               ("FM", FM, "i"),
                               ("BS", descDict1["BS"], "Z"),
                               ("FP", int("Pass" in descDict1["FP"]), "i"),
                               ("PV", descDict1["PV"], "Z"),
                               ("FA", descDict1["FA"], "Z"),
                               ("ND", ND2, "i"),
                               ("NF", odiv(ND2, float(FM)), "f"),
                               ("RG", "default", "Z"),
                               ("AF", r2FracAlign, "f"),
                               ("SF", r2FracSC, "f")])
        else:
            read1bam.set_tags([("RP", coorString, "Z"),
                               ("SC", contigSetStr, "Z"),
                               ("FM", FM, "i"),
                               ("BS", descDict1["BS"], "Z"),
                               ("FP", int("Pass" in descDict1["FP"]), "i"),
                               ("PV", descDict1["PV"], "Z"),
                               ("FA", descDict1["FA"], "Z"),
                               ("ND", ND1, "i"),
                               ("NF", 1. * ND1 / FM, "f"),
                               ("AF", r1FracAlign, "f"),
                               ("SF", r1FracSC, "f")])
            read2bam.set_tags([("RP", coorString, "Z"),
                               ("SC", contigSetStr, "Z"),
                               ("FM", FM, "i"),
                               ("BS", descDict1["BS"], "Z"),
                               ("FP", int("Pass" in descDict1["FP"]), "i"),
                               ("PV", descDict1["PV"], "Z"),
                               ("FA", descDict1["FA"], "Z"),
                               ("ND", ND2, "i"),
                               ("NF", 1. * ND2 / FM, "f"),
                               ("AF", r2FracAlign, "f"),
                               ("SF", r2FracSC, "f")])
        # I used to mark the BAMs at this stage, but it's not appropriate to
        # do so until after indel realignment.
        obw(read1bam)
        obw(read2bam)
    suppBAM.close()
    outBAM.close()
    postFilterBAM.close()
    return outBAMFile


def compareRecs(RecordList):
    Success = True
    seqs = list(cmap(oag("seq"), RecordList))
    seqs = [str(record.seq) for record in RecordList]
    stackArrays = tuple([npchararray(s, itemsize=1) for s in seqs])
    seqArray = nvstack(stackArrays)
    # print(repr(seqArray))

    quals = nparray(list(cmap(oag("query_qualities"),
                         RecordList)))
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
    newSeq = "".join([letterNumDict[i] for i in nargmax(qualAllSum, 0)])
    MaxPhredSum = np.amax(qualAllSum, 0)  # Avoid calculating twice.
    phredQuals = nsubtract(nmultiply(2, MaxPhredSum),
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
    inputHandle = pysam.Samfile(inBAM, 'rb')
    outputHandle = pysam.Samfile(outBAM, 'wb', template=inputHandle)
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


def criteriaTest(read1, read2, filterSet="default", minFamSize=3):
    """
    Tool for filtering a pair of BAM files by criteria.
    Note: complexity filter is not needed for shades protocol.
    """
    list = ("adapter barcode complexity editdistance family "
            "ismapped qc notinbed").split()
    Logger = logging.getLogger("Primarylogger")
    try:
        assert isinstance(filterSet, str)  # This should be a string.
    except AssertionError:
        if(isinstance(filterSet, list)):
            pl("filterSet is a list, which was not expected. repr: {}".format(
                repr(filterSet)))
            pl("Joining filterSet into a string.")
            filterSet = ''.join(filterSet)
            if(isinstance(filterSet, str) is False):
                HTSUtils.FacePalm("filterSet isn't a string after joining.")
        else:
            pl("Filter set is not as expected! Repr: {}".format(
                repr(filterSet)))
            raise ValueError("What is this? filterSet is not string or list.")
    if(filterSet == "default"):
        pl("List of valid filters: {}".format(', '.join(list)))
        raise ValueError("Filter must be set!")
    filterSet = filterSet.split(',')
    for filt in filterSet:
        if filt not in list:
            pl("filterSet provided: {}".format(filterSet))
            raise ValueError("Select valid filter(s) - {}".format(list))

    if("adapter" in filterSet):
        ALValue1 = int(read1.opt("FP"))
        ALValue2 = int(read2.opt("FP"))
        if(sum([ALValue1, ALValue2]) != 2):
            return False

    if("barcode" in filterSet):
        if(read1.opt("BS") != read2.opt("BS")):
            Logger.debug(
                "Barcode sequence didn't match. Are you running shades?")
            return False

    if("editdistance" in filterSet):
        NMValue1 = int(read1.opt("NM"))
        NMValue2 = int(read2.opt("NM"))
        if(NMValue1 == 0 and NMValue2 == 0):
            return False

    if("family" in filterSet):
        FMValue1 = int(read1.opt("FM"))
        FMValue2 = int(read2.opt("FM"))
        if(FMValue1 < int(minFamSize) or FMValue2 < int(minFamSize)):
            Logger.debug(("This family didn't survive. Their FMValues:" +
                          "{}, {}".format(FMValue1, FMValue2)))
            return False

    if("ismapped" in filterSet):
        if(read1.is_unmapped and read2.is_unmapped):
            return False

    if("qc" in filterSet):
        if(read1.is_qcfail or read2.is_qcfail):
            return False

    return True


def GenerateFamilyHistochart(BCIdx, output="default"):
    if(output == "default"):
        output = '.'.join(BCIdx.split('.')[:-1]) + '.hist.txt'
    Str = ("cat {} | awk '{{print $1}}' | sort | uniq -c | ".format(BCIdx) +
           "awk 'BEGIN {{OFS=\"\t\"}};{{print $1,$2}}' | sort " +
           "-k1,1n > {}".format(output))
    pl("Command str: {}".format(Str))
    HTSUtils.PipedShellCall(Str)
    pl("Family size histochart: {}".format(output))
    return output


def pairedFilterBam(inputBAM, passBAM="default",
                    failBAM="default", criteria="default",
                    deleteFailBam=False):
    """
    Filters out both reads in a pair based on a list of ","-separated criteria.
    Both reads must pass for the pair to be written
    Required: [sb]am file, name-sorted, keep unmapped reads,
    and remove supplementary and secondary alignments.
    criteria must be a string. If multiple filters are set,
    they must be comma-separated and have no spaces.
    """
    cStr = ("pairedFilterBam({}, passBAM='{}', ".format(inputBAM, passBAM) +
            "failBAM='{}', criteria='{}', delete".format(failBAM, criteria) +
            "FailBam='{}')".format(deleteFailBam))
    pl("Command required to reproduce this call: {}".format(cStr))
    if(criteria == "default"):
        raise NameError("Filter Failed: Criterion Not Set.")
    if(passBAM == "default"):
        passBAM = '.'.join(inputBAM.split('.')[0:-1])
        passBAM += ".{}P.bam".format(criteria)
    if(failBAM == "default"):
        failBAM = '.'.join(inputBAM.split('.')[0:-1])
        failBAM += ".{}F.bam".format(criteria)
    pl("pairedFilterBam. Input: {}. Pass: {}".format(inputBAM, passBAM))
    inBAM = pysam.Samfile(inputBAM, "rb")
    passFilter = pysam.Samfile(passBAM, "wb", template=inBAM)
    failFilter = pysam.Samfile(failBAM, "wb", template=inBAM)
    criteriaList = criteria.lower()
    pl("Criteria string is: {}".format(criteriaList))
    pl("Following criteria will be used: {}".format(
        str(criteriaList.split(','))))
    for read in inBAM:
        failed = False
        if(read.is_read1):
            read1 = read
            continue
        if(read.is_read2):
            read2 = read
            try:
                assert read1.qname == read2.qname  # Sanity check
            except AssertionError:
                pl("Failed sanity check. Is the alignment file name-sorted?")
                raise ThisIsMadness("Please inspect your bam file!")
            result = criteriaTest(read1, read2, filterSet=criteriaList)
            if(result is False):
                failed = True
                failFilter.write(read1)
                failFilter.write(read2)
                continue
            if(failed):
                continue
            passFilter.write(read1)
            passFilter.write(read2)
    passFilter.close()
    failFilter.close()
    inBAM.close()
    if(str(deleteFailBam).lower() == "true"):
        os.remove(failBAM)
    return passBAM, failBAM


def removeSecondary(inBAM, outBAM="default"):
    if(outBAM == "default"):
        outBAM = '.'.join(inBAM.split('.')[0:-1]) + '.2ndrm.bam'
    input = pysam.Samfile(inBAM, "rb")
    output = pysam.Samfile(outBAM, "wb", template=input)
    pl(("Attempting to remove secondary"))
    for entry in input:
        if(entry.is_secondary or entry.flag >> 11 == 1):
            continue
        else:
            output.write(entry)
    return outBAM


def Sam2Bam(insam, outBAM):
    pl("Sam2Bam converting {} to {}".format(insam, outBAM))
    output = open(outBAM, 'w', 0)
    command_str = 'samtools view -Sbh {}'.format(insam)
    pl((command_str))
    subprocess.check_call(shlex.split(command_str), stdout=output, shell=False)
    return(command_str, outBAM)


# Taken from Scrutils, by Jacob Durtschi
def SamtoolsBam2fq(bamPath, outFastqPath):
    pl("SamtoolsBam2fq converting {}".format(bamPath))
    # Build commands that will be piped
    samtoolsBam2fqCommand = ['samtools', 'bam2fq', bamPath]
    gzipCommand = ['gzip']

    # Open output fastq.gz file to be piped into
    outFastq = open(outFastqPath, 'w')

    # Call piped commands
    process1 = subprocess.Popen(samtoolsBam2fqCommand,
                                stdout=subprocess.PIPE, shell=False)
    process2 = subprocess.Popen(gzipCommand, stdin=process1.stdout,
                                stdout=outFastq, shell=False)
    # process1.stdout.close()
    # process1.wait()
    process2.wait()
    outFastq.flush()
    outFastq.close()
    # if processErr:

    return outFastqPath


def singleBarcodeTagging(fastq, bam, outputBAM="default", suppBam="default"):
    if(outputBAM == "default"):
        outputBAM = '.'.join(bam.split('.')[0:-1]) + ".tagged.bam"
    if(suppBam == "default"):
        suppBam = bam.split('.')[0] + '.2ndSupp.bam'
    pl("singleBarcodeTagging. Fq: {}. outputBAM: {}".format(bam, outputBAM))
    reads = SeqIO.parse(fastq, "fastq")
    # inBAM = removeSecondary(args.bam_file) #Artefactual code
    postFilterBAM = pysam.Samfile(bam, "rb")
    suppBAM = pysam.Samfile(suppBam, "wb", template=postFilterBAM)
    outBAM = pysam.Samfile(outputBAM, "wb", template=postFilterBAM)
    for entry in postFilterBAM:
        if(entry.is_secondary or entry.flag >> 11 == 1):
            suppBAM.write(entry)
            continue
        else:
            try:
                tempRead = reads.next()
            except StopIteration:
                break
        descDict = getdesc(tempRead.description)
        for key in descDict.keys():
            entry.setTag(key, descDict[key])
        if("Pass" not in descDict["FP"]):
            pl(("Standard filter for barcode failed! Val: {}".format(
                descDict["FP"])))
            pl(("Tags dictionary is {}".format(descDict)))
            raise ValueError("Something has gone wrong!")
        if("Pass" in descDict["FP"]):
            entry.tags = entry.tags + [("FP", 1)]
        else:
            entry.tags = entry.tags + [("FP", 0)]
        outBAM.write(entry)
    outBAM.close()
    postFilterBAM.close()
    return outputBAM


def SingleConsolidate(inBAM, outBAM="default"):
    if(outBAM == "default"):
        outBAM = '.'.join(inBAM.split('.')[0:-1]) + "consolidated.bam"
    pl("SingleConsolidate. Input: {}. Output: {}".format(inBAM, outBAM))
    inputHandle = pysam.Samfile(inBAM, 'rb')
    outputHandle = pysam.Samfile(outBAM, 'wb', template=inputHandle)
    workBC = ""
    Set = []
    for record in inputHandle:
        barcodeRecord = [
            tagSet for tagSet in record.tags if tagSet[0] == "BS"][0][1]
        # print("Working: {}. Current: {}.".format(workBC,barcodeRecord))
        # print("name of read with this barcode: {}".format(record.qname))
        # print("Working set: {}".format(Set))
        if(workBC == ""):
            workBC = barcodeRecord
            Set = []
            Set.append(record)
            continue
        elif(workBC == barcodeRecord):
            # print(record.qname)
            # print(workBC==barcodeRecord)
            # print("WorkBC: {}. Current: {}.".format(workBC,barcodeRecord))
            try:
                assert([tagSet for tagSet in Set[0].tags if tagSet[
                    0] == "BS"][0][1] == [
                    tagSet for tagSet in record.tags if tagSet[
                        0] == "BS"][0][1])
            except AssertionError:
                BS1 = [tagSet for tagSet in Set[
                    0].tags if tagSet[0] == "BS"][0][1]
                BS2 = [tagSet for tagSet in record.tags if tagSet[
                    0] == "BS"][0][1]
                pl("BS1: {}. BS2: {}".format(BS1, BS2))
                raise AssertionError("Well, there you go.")
            Set.append(record)
            continue
        elif(workBC != barcodeRecord):
            mergeRec, success = compareRecs(Set)
            if(success):
                outputHandle.write(mergeRec)
            Set = []
            Set.append(record)
            workBC = barcodeRecord
            continue
        else:
            raise RuntimeError("This code should be unreachable")
    inputHandle.close()
    outputHandle.close()
    return outBAM


@cython.returns(cython.bint)
def singleCriteriaTest(read, filter="default"):
    list = "adapter complexity editdistance family ismapped qc".split(' ')

    if(filter == "default"):
        pl(("List of valid filters: {}".format(', '.join(list))))
        raise ValueError("Filter must be set!.")

    if(filter not in list):
        raise ValueError("Filter not supported. The list is: {}".format(list))

    if(filter == "adapter"):
        ALloc1 = [i for i, j in enumerate(read.tags) if j[0] == "FP"][0]

        try:
            ALValue1 = int(read.tags[ALloc1][1])
        except IndexError:
            ALValue1 = 0
        if(ALValue1 == 0):
            return False

    if(filter == "complexity"):
        a = "AAAAAAAAAA"
        t = "TTTTTTTTTT"
        c = "CCCCCCCCCC"
        g = "GGGGGGGGGG"
        rt = str(read.tags)
        if(a in rt or t in rt or g in rt or c in rt):
            return False

    if(filter == "editdistance"):
        NMloc1 = [i for i, j in enumerate(read.tags) if j[0] == "NM"][0]
        if(read.tags[NMloc1][1] == 0):
            return False

    if(filter == "family"):
        FMloc1 = [i for i, j in enumerate(read.tags) if j[0] == "FM"][0]
        FMValue1 = int(read.tags[FMloc1][1])
        if(FMValue1 < 3):
            return False

    if(filter == "ismapped"):
        if(read.is_unmapped):
            return False

    if(filter == "qc"):
        if(read.is_qcfail):
            return False

    return True


def singleFilterBam(inputBAM, passBAM="default",
                    failBAM="default", criteria="default"):
    """
    Filters out both reads in a pair based on a list of
    comma-separated criteria. Both reads must pass to be written.
    Required: [SB]AM file, coordinate-sorted,
    supplementary and secondary alignments removed,unmapped reads retained.
    """
    cdef pysam.calignmentfile.AlignedSegment read
    if(criteria == "default"):
        raise NameError("Filter Failed: Criterion Not Set.")
    if(passBAM == "default"):
        passBAM = '.'.join(inputBAM.split('.')[0:-1])
        passBAM += ".{}P.bam".format(criteria)
    if(failBAM == "default"):
        failBAM = '.'.join(inputBAM.split('.')[0:-1])
        failBAM += ".{}F.bam".format(criteria)
    pl("singleFilterBam. Input: {}. Pass: {}".format(inputBAM, passBAM))
    inBAM = pysam.Samfile(inputBAM, "rb")
    passFilter = pysam.Samfile(passBAM, "wb", template=inBAM)
    failFilter = pysam.Samfile(failBAM, "wb", template=inBAM)
    criteriaList = criteria.lower().split(',')
    for i, entry in enumerate(criteriaList):
        pl(("Criteria #{} is \"{}\"".format(i, entry)))
    for read in inBAM:
        failed = False
        for criterion in criteriaList:
            result = singleCriteriaTest(read, filter=criterion)
            if(result is False):
                failFilter.write(read)
                failed = True
                break
            else:
                continue
        if(failed):
            continue
        passFilter.write(read)
    passFilter.close()
    failFilter.close()
    inBAM.close()
    return passBAM, failBAM


def GetRPsWithI(inBAM, outBAM="default"):
    """
    Gets read pairs with an I in both R1 and R2's cigar strings.
    Must be namesorted!
    If outBAM is left as default, it chooses to write to a filename based
    on the inBAM name.
    """
    cdef pysam.calignmentfile.AlignedSegment read1
    cdef pysam.calignmentfile.AlignedSegment read2
    cdef pysam.calignmentfile.AlignedSegment entry
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    if(outBAM == "default"):
        outBAM = ".".join(inBAM.split(".")[0:-1]) + ".InsertedReadPairs.bam"
    outHandle = pysam.AlignmentFile(outBAM, "wb", template=inHandle)
    ohw = outHandle.write
    for entry in inHandle:
        if(entry.is_read1):
            read1 = entry
            continue
        if(entry.is_read2):
            read2 = entry
        try:
            assert read1.query_name == read2.query_name
        except AssertionError:
            HTSUtils.FacePalm("Input fastq is not name sorted or is off in "
                              "some other way! Abort!")
        if(read2.cigarstring is None or read1.cigarstring is None or
           "I" not in read1.cigarstring or "I" not in read2.cigarstring):
            continue
        ohw(read1)
        ohw(read2)
    outHandle.close()
    return outBAM
