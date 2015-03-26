import re
import subprocess
import time
import logging

import cython
import numpy as np

from MawCluster import BCBam, VCFWriters, BCVCF
from MawCluster import BCFastq
from utilBMF import HTSUtils
from MawCluster import PileupUtils
from MawCluster.SVUtils import GetSVRelevantRecordsPaired as SVRP
from utilBMF.HTSUtils import printlog as pl, ThisIsMadness
from MawCluster.BCVCF import VCFStats
from MawCluster.FFPE import GetDeaminationFrequencies, FilterByDeaminationFreq


@cython.locals(calcCoverage=cython.bint, coverageForAllRegions=cython.bint,
               addRG=cython.bint, mincov=cython.long, rLen=cython.long)
def pairedBamProc(consfq1, consfq2, consfqSingle="default", opts="",
                  bamPrefix="default", ref="default", aligner="default",
                  barIndex="default",
                  bed="/yggdrasil/workspace/Barcode_Noah/cfDNA_targets.bed",
                  mincov=5,
                  abrapath="default",
                  coverageForAllRegions=False,
                  calcCoverage=True,
                  bwapath="default",
                  picardPath="default",
                  addRG=False, PL="ILLUMINA",
                  SM="default", CN="default",
                  RG="default", ID="default",
                  realigner="gatk", gatkpath="default", dbsnp="default",
                  rLen=-1, intelDeflator="default"):
    """
    Performs alignment and sam tagging of consolidated fastq files.
    Note: the i5/i7 indexing strategy ("Shades") does not use the consfqSingle
    """
    if(ref == "default"):
        raise ValueError("Reference index required!")
    if(bamPrefix == "default"):
        bamPrefix = '.'.join(consfq1.split('.')[0:-1])
    if(aligner == "default"):
        pl("No aligner set, defaulting to bwa mem.")
        aligner = "mem"
    if(aligner == "mem"):
        if(addRG is False):
            outBAMProperPair = HTSUtils.align_bwa_mem(
                consfq1, consfq2, ref=ref, opts=opts, path=bwapath)
        else:
            outBAMProperPair = HTSUtils.align_bwa_mem_addRG(
                consfq1, consfq2, ref=ref, opts=opts, path=bwapath,
                RG=RG, ID=ID, CN=CN, PL=PL, SM=SM)
        if(consfqSingle != "default"):
            HTSUtils.FacePalm("This step is not relevant to shades.")
    elif(aligner == "aln"):
        if(addRG is False):
            outBAMProperPair = HTSUtils.align_bwa_aln(consfq1, consfq2,
                                                      ref=ref, opts=opts)
        else:
            outBAMProperPair = HTSUtils.align_bwa_aln_addRG(
                consfq1, consfq2, ref=ref, opts=opts, RG=RG, SM=SM,
                CN=CN, PL=PL, picardPath=picardPath, ID=ID)
        if(consfqSingle != "default"):
            HTSUtils.FacePalm("This step is not required "
                              "or important for shades.")
    else:
        raise ValueError("Sorry, only bwa is supported currently.")
    if(rLen < 0):
        pl("rLen < 0 in pairedBamproc. This typically means that this "
           "is not set.")
    if(picardPath == "default"):
        pl("Warning: path to picard jar not set. This isn't required for much"
           ", but in case something dies later, this could be responsible",
           level=logging.DEBUG)
    pl("Now tagging BAM with custom SAM tags.")
    taggedBAM = BCBam.pairedBarcodeTagging(
        consfq1, consfq2, outBAMProperPair, bedfile=bed)
    pl("Now splitting the BAM into read 1 and read 2 files.")
    pl("Now generating double barcode index.")
    if(realigner == "abra"):
        if(abrapath == "default"):
            raise ThisIsMadness("abrapath must be set for abra to be the "
                                "realigner")
        coorSortFull = HTSUtils.CoorSortAndIndexBam(taggedBAM)
        realignedFull = BCBam.AbraCadabra(coorSortFull, ref=ref, bed=bed,
                                          jar=abrapath, rLen=rLen,
                                          intelPath=intelDeflator)
    elif(realigner == "gatk"):
        coorSortFull = HTSUtils.CoorSortAndIndexBam(taggedBAM)
        realignedFull = BCBam.GATKIndelRealignment(coorSortFull, ref=ref,
                                                   bed=bed, gatk=gatkpath,
                                                   dbsnp=dbsnp)
    else:
        pl("ABRA path not provided. Skipping realignment.")
        realignedFull = taggedBAM
        # realignedFull = BCBam.AbraCadabra(taggedBAM, ref=ref, bed=bed)
    namesortedRealignedFull = HTSUtils.NameSort(realignedFull, uuid=True)
    if(barIndex != "default"):
        p = subprocess.Popen(["wc", "-l", barIndex], stdout=subprocess.PIPE)
        out, err = p.communicate()
        pl("Number of families found: {}".format(
            re.findall(r'\d+', out)[0]))
        histochart = BCBam.GenerateFamilyHistochart(barIndex)
        pl("Histochart of family sizes: {}".format(histochart))
    # UNCOMMENT THIS BLOCK IF YOU WANT TO START MESSING WITH RESCUE
    """
        pl("Rescue step, marking the BD as their Hamming distance.")
        newRef = GenerateBarcodeIndexReference(uniqueBigFamilies)
        indexBowtie(newRef)
        mergedFastq = mergeSequencesFastq(tags1, tags2,)
        joiningSAM = CustomRefBowtiePaired(mergedFastq,newRef)
        return
        joinedFamilies = fuzzyJoining(familyMarked,joiningSAM)
        pl("joinedFamilies is {}".format(joinedFamilies))
    """
    # This step not needed for shades protocol, as fastq
    # families have already been filtered for size.
    # familyP, familyF = BCBam.pairedFilterBam(
    #    families, criteria="family")
    tempBAMPrefix = '.'.join(namesortedRealignedFull.split('.')[0:-1])
    summary = ".".join(namesortedRealignedFull.split('.') + ['SV', 'txt'])

    SVBam, MarkedFamilies = SVRP(namesortedRealignedFull,
                                 bedfile=bed,
                                 tempBAMPrefix=tempBAMPrefix,
                                 summary=summary)
    if(SVBam != "NotWritten"):
        pl(("{} is the bam with all reads considered relevant ".format(SVBam) +
            "to structural variants."))
    else:
        pl("Structural variant bam not written. Ask for it next time!")
    pl("All records BAM: %s" % MarkedFamilies)
    # SVOutputFile = BCBam.CallTranslocations(SVBam, bedfile=bed)
    pl("Change of plans - now, the SV-marked BAM is not used for "
       "SNP calling due to the differing alignment needs.")
    coorSorted = HTSUtils.CoorSortAndIndexBam(MarkedFamilies)
    return coorSorted


@cython.locals(overlapLen=cython.long,
               stringency=cython.float)
def pairedFastqShades(inFastq1, inFastq2, indexfq="default", stringency=0.9,
                      p3Seq="default", p5Seq="default",
                      overlapLen=6):
    pl("Beginning pairedFastqShades for {}, {}".format(inFastq1, inFastq2))
    bcFastq1, bcFastq2 = BCFastq.FastqPairedShading(inFastq1,
                                                    inFastq2,
                                                    indexfq=indexfq)
    BSortFq1, BSortFq2 = BCFastq.BarcodeSortBoth(bcFastq1, bcFastq2)
    BConsFastq1, BConsFastq2 = BCFastq.pairedFastqConsolidateFaster(
        BSortFq1, BSortFq2, stringency=0.9)
    if(p3Seq != "default"):
        BConsFastq1, BConsFastq2 = BCFastq.CutadaptPaired(
            BConsFastq1, BConsFastq2, overlapLen=overlapLen,
            p3Seq=p3Seq, p5Seq=p5Seq)
    return BConsFastq1, BConsFastq2


@cython.locals(minMQ=cython.long, minBQ=cython.long, MakeVCF=cython.bint,
               MakeCoverageBed=cython.bint, MakePileupTsv=cython.bint,
               minFA=cython.long, minFracAgreed=cython.float,
               deaminationPVal=cython.float)
def pairedVCFProc(consMergeSortBAM,
                  ref="default",
                  opts="",
                  bed="default",
                  minMQ=10,
                  minBQ=20,
                  MakePileupTsv=False,
                  MakeVCF=True,
                  MakeCoverageBed=False,
                  commandStr="default",
                  minFA=2, minFracAgreed=0.667,
                  exp="", deaminationPVal=0.05):
    """
    Lumps together VCF processing.
    exp is a string from a comma-joined list of strings.
    If "ffpe" is in exp.lower(), then an FFPE deamination step will
    be run to filter out those artefacts.
    """
    if(bed == "default"):
        raise ValueError("Bed file location must be set!")
    if(ref == "default"):
        raise ValueError("Reference index location must be set!")
    # Consolidating families into single reads
    # Variant Calling Step using MPileup
    # print("Now filtering for reads with NM > 0 only if you want to.")
    Results = {}
    if(MakeCoverageBed):
        OutBed = PileupUtils.CalcWithinBedCoverage(consMergeSortBAM,
                                                   bed=bed,
                                                   minMQ=minMQ,
                                                   minBQ=minBQ)
        Results["bed"] = OutBed
    if(MakePileupTsv):
        PileupTSV = PileupUtils.CustomPileupToTsv(consMergeSortBAM,
                                                  bedfile=bed,
                                                  minMQ=minMQ,
                                                  minBQ=minBQ)
        pl("PileupTSV: {}".format(PileupTSV))
        Results["tsv"] = PileupTSV
    if(MakeVCF):
        SNP_VCF = VCFWriters.SNVCrawler(consMergeSortBAM,
                                        minMQ=minMQ,
                                        minBQ=minBQ,
                                        reference=ref,
                                        commandStr=commandStr,
                                        reference_is_path=True,
                                        bed=bed, minFA=minFA,
                                        minFracAgreed=minFracAgreed,
                                        experiment=exp)
        CleanedVCF = BCVCF.FilterVCFFileByBed(SNP_VCF, bed)
        if("ffpe" in exp.lower()):
            deaminationFrequency = np.mean(
                GetDeaminationFrequencies(CleanedVCF), dtype=np.longdouble)
            FFPEFilteredVCF = FilterByDeaminationFreq(
                CleanedVCF, ctfreq=deaminationFrequency,
                pVal=deaminationPVal)
            finalVCF = FFPEFilteredVCF
        else:
            finalVCF = CleanedVCF
        pl("SNP VCF: {}".format(finalVCF))
        Results["vcf"] = finalVCF
        VCFStatsFile = VCFStats(finalVCF)
        Results["vcfstats"] = VCFStatsFile
    # AlleleFreqTSV = PileupUtils.AlleleFrequenciesByBase(consMergeSortBAM,
    #                                                     bedfile=bed)
    # This is probably useless given that I'm doing this "manually",
    # but I'm keeping this in here for good measure.
    return Results
