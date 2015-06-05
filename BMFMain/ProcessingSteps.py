try:
    import re2 as re
except ImportError:
    import re
import subprocess
import time
import logging
from subprocess import check_call, CalledProcessError

import cython
import numpy as np

from MawCluster import BCBam, VCFWriters, BCVCF
from MawCluster import BCFastq
from MawCluster.BCFastq import TrimHomingPaired
from utilBMF import HTSUtils
from MawCluster import PileupUtils
from MawCluster.SVUtils import GetSVRelevantRecordsPaired as SVRP
from utilBMF.HTSUtils import printlog as pl, TrimExt
from utilBMF.QC import GetAllQCMetrics, GetFamSizeStats
from utilBMF.ErrorHandling import ThisIsMadness as Tim
from MawCluster.BCVCF import VCFStats
from MawCluster.FFPE import GetDeaminationFrequencies, FilterByDeaminationFreq


@cython.locals(calcCoverage=cython.bint, coverageForAllRegions=cython.bint,
               addRG=cython.bint, mincov=int, rLen=int)
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
                  realigner="abra", gatkpath="default", dbsnp="default",
                  rLen=-1, intelDeflator="default", minAF=0.0):
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
    if("gatk" in realigner.lower()):
        pl("Since realigning with GATK IndelRealigner, "
           "addRG is being set to True.")
        addRG = True
    if(aligner == "mem"):
        outBAMProperPair = BCBam.AlignAndTagMem(
            consfq1, consfq2, ref=ref, opts=opts, addRG=addRG)
    elif(aligner == "aln"):
        if(addRG is False):
            outBAMProperPair = HTSUtils.align_bwa_aln(consfq1, consfq2,
                                                      ref=ref, opts=opts)
        else:
            outBAMProperPair = HTSUtils.align_bwa_aln_addRG(
                consfq1, consfq2, ref=ref, opts=opts, RG=RG, SM=SM,
                CN=CN, PL=PL, picardPath=picardPath, ID=ID)
        if(consfqSingle != "default"):
            raise Tim("This step is not required or important for shades.")
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
        consfq1, consfq2, outBAMProperPair, bedfile=bed, realigner=realigner,
        ref=ref, minAF=minAF)
    # check_call(["rm", outBAMProperPair])
    pl("Now realigning with: %s" % realigner)
    if("abra" in realigner.lower()):
        if(abrapath == "default"):
            raise Tim("abrapath must be set for abra to be the realigner")
        coorSortFull = HTSUtils.CoorSortAndIndexBam(taggedBAM)
        realignedFull = BCBam.AbraCadabra(coorSortFull, ref=ref, bed=bed,
                                          jar=abrapath, rLen=rLen,
                                          intelPath=intelDeflator)
        if("gatk" in realigner.lower()):
            pl("Realigning around known indels, too. "
               "Lots of steps, hard to know which ones matter.")
            coorSortFull = HTSUtils.CoorSortAndIndexBam(realignedFull)
            realignedFull = BCBam.GATKIndelRealignment(coorSortFull, ref=ref,
                                                       bed=bed, gatk=gatkpath,
                                                       dbsnp=dbsnp)
    elif("gatk" in realigner.lower().split(",")):
        coorSortFull = HTSUtils.CoorSortAndIndexBam(taggedBAM)
        realignedFull = BCBam.GATKIndelRealignment(coorSortFull, ref=ref,
                                                   bed=bed, gatk=gatkpath,
                                                   dbsnp=dbsnp)
    else:
        pl("Skipping realignment.")
        realignedFull = taggedBAM
        # realignedFull = BCBam.AbraCadabra(taggedBAM, ref=ref, bed=bed)
    namesortedRealignedFull = HTSUtils.NameSort(realignedFull, uuid=True)
    # check_call(["rm", realignedFull])
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
    check_call(["rm", namesortedRealignedFull])
    coorSorted = HTSUtils.CoorSortAndIndexBam(MarkedFamilies)
    check_call(["rm", MarkedFamilies])
    return coorSorted


@cython.locals(overlapLen=int,
               stringency=float)
def pairedFastqShades(inFastq1, inFastq2, indexfq="default", stringency=0.95,
                      p3Seq="default", p5Seq="default",
                      overlapLen=6, sortMem="6G", inline_barcodes=False,
                      homing=None, bcLen=-1, head=0):
    pl("Beginning pairedFastqShades for {}, {}".format(inFastq1, inFastq2))
    if(inline_barcodes is False):
        bcFastq1, bcFastq2 = BCFastq.FastqPairedShading(inFastq1,
                                                        inFastq2,
                                                        indexfq=indexfq,
                                                        head=head)
    else:
        bcFastq1, bcFastq2 = TrimHomingPaired(inFastq1, inFastq2,
                                              homing=homing, bcLen=bcLen)
    BSortFq1, BSortFq2 = BCFastq.BarcodeSortBoth(bcFastq1, bcFastq2,
                                                 sortMem=sortMem)
    # check_call(["rm", bcFastq1, bcFastq2])
    BConsFastq1, BConsFastq2 = BCFastq.pairedFastqConsolidate(
        BSortFq1, BSortFq2, stringency=0.9)
    pl("Parameters for cutadapt are p3Seq={}, p5Seq={}".format(p3Seq, p5Seq))
    if(p3Seq != "default"):
        pl("Running cutadapt ...")
        BConsFastq1, BConsFastq2 = BCFastq.CutadaptPaired(
            BConsFastq1, BConsFastq2, overlapLen=overlapLen,
            p3Seq=p3Seq, p5Seq=p5Seq)
    else:
        pl("Skipping cutadapt ...")
    check_call(["rm", BSortFq1, BSortFq2])
    famStats = GetFamSizeStats(
        BConsFastq1,
        outfile=TrimExt(inFastq1) + ".famstats.txt")
    return BConsFastq1, BConsFastq2


@cython.locals(minMQ=int, minBQ=int, MakeVCF=cython.bint,
               MakeCoverageBed=cython.bint, MakePileupTsv=cython.bint,
               minFA=int, minFracAgreed=float,
               deaminationPVal=float)
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
                  exp="", deaminationPVal=0.05,
                  conf="default", parallel=True):
    """
    Lumps together VCF processing.
    exp is a string from a comma-joined list of strings.
    If "ffpe" is in exp.lower(), then an FFPE deamination step will
    be run to filter out those artefacts.
    """
    if(bed == "default"):
        raise ValueError("Bed file location must be set!")
    if(conf == "default"):
        raise ValueError("config file location must be set!")
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
        if(parallel):
            CleanedVCF = VCFWriters.PSNVCall(consMergeSortBAM,
                                             conf=conf)
        else:
            CleanedVCF = VCFWriters.SNVCrawler(consMergeSortBAM, conf=conf)
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
