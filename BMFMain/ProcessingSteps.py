import re
import subprocess

from MawCluster import BCBam, VCFWriters, BCVCF
from MawCluster import BCFastq
from utilBMF import HTSUtils
from MawCluster import PileupUtils
from MawCluster.SVUtils import GetSVRelevantRecordsPaired as SVRP
from utilBMF.HTSUtils import printlog as pl
from MawCluster.BCVCF import VCFStats


def pairedBamProc(consfq1, consfq2, consfqSingle="default", opts="",
                  bamPrefix="default", ref="default", aligner="default",
                  barIndex="default",
                  bed="/yggdrasil/workspace/Barcode_Noah/cfDNA_targets.bed",
                  mincov=5,
                  abrapath="default",
                  coverageForAllRegions=False):
    """
    Performs alignment and sam tagging of consolidated fastq files.
    Note: the i5/i7 indexing strategy ("Shades") does not use the consfqSingle
    """
    if(ref == "default"):
        raise ValueError("Reference index required!")
    if(barIndex == "default"):
        raise ValueError(("Barcode index required - generate one at the "
                          "pre-merge fastq stage."))
    if(bamPrefix == "default"):
        bamPrefix = '.'.join(consfq1.split('.')[0:-1])
    if(aligner == "default"):
        pl("No aligner set, defaulting to bwa mem.")
        aligner = "mem"
    if(aligner == "mem"):
        outBAMProperPair = HTSUtils.align_bwa_mem(
            consfq1, consfq2, ref=ref, opts=opts)
        if(consfqSingle != "default"):
            HTSUtils.FacePalm("This step is not required "
                              "or important for shades.")
    elif(aligner == "aln"):
        outBAMProperPair = HTSUtils.align_bwa_aln(consfq1, consfq2, ref=ref,
                                                  opts=opts)
        if(consfqSingle != "default"):
            HTSUtils.FacePalm("This step is not required "
                              "or important for shades.")
    else:
        raise ValueError("Sorry, only bwa is supported currently.")
    pl("Now tagging BAM with custom SAM tags.")
    taggedBAM = BCBam.pairedBarcodeTagging(
        consfq1, consfq2, outBAMProperPair)
    pl("Now splitting the BAM into read 1 and read 2 files.")
    pl("Now generating double barcode index.")
    if(abrapath != "default"):
        realignedFull = BCBam.AbraCadabra(taggedBAM, ref=ref, bed=bed,
                                          jar=abrapath)
    else:
        pl("ABRA path not provided. Skipping realignment.")
        realignedFull = taggedBAM
        # realignedFull = BCBam.AbraCadabra(taggedBAM, ref=ref, bed=bed)
    namesortedRealignedFull = HTSUtils.NameSort(realignedFull, uuid=True)
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
    SVBam, MarkedFamilies = SVRP(namesortedRealignedFull,
                                 bedfile=bed,
                                 tempBAMPrefix=namesortedRealignedFull[0:-4],
                                 summary=(namesortedRealignedFull[0:-4] +
                                          '.SV.txt'))
    pl(("{} is the bam with all reads considered relevant ".format(SVBam) +
        "to structural variants."))
    # SVOutputFile = BCBam.CallTranslocations(SVBam, bedfile=bed)
    pl("Change of plans - now, the SV-marked BAM is not used for "
       "SNP calling due to an error.")
    coorSorted = BCBam.CoorSort(namesortedRealignedFull)
    if(coverageForAllRegions is False):
        CoverageBed = PileupUtils.CalcWithinBedCoverage(coorSorted, bed=bed)
    else:
        CoverageBed = PileupUtils.BamToCoverageBed(coorSorted, mincov=mincov)
    pl("Coverage bed: {}".format(CoverageBed))
    return coorSorted


def pairedFastqShades(inFastq1, inFastq2, indexfq="default", stringency=0.75):
    bcFastq1, bcFastq2 = BCFastq.FastqPairedShading(inFastq1,
                                                    inFastq2,
                                                    indexfq=indexfq,
                                                    gzip=False)
    if(indexfq == "default"):
        HTSUtils.FacePalm("pairedFastqShades requires an index fastq.")
    pl("Beginning pairedFastqShades for {}, {}".format(inFastq1, inFastq2))
    barcodeIndex = BCFastq.GenerateShadesIndex(indexfq)
    (FamFqs, SingleFqs, numReads,
     numReadsWFam) = BCFastq.GetFamilySizePaired(bcFastq1,
                                                 bcFastq2, barcodeIndex)
    FamFq1 = FamFqs[0]
    FamFq2 = FamFqs[1]
    SingleFq1 = SingleFqs[0]
    SingleFq2 = SingleFqs[1]
    pl("Number of reads total: " + str(numReads))
    pl("Number of reads with >=3 family members: " + str(numReadsWFam))
    BSortFq1 = BCFastq.BarcodeSort(FamFq1)
    BSortFq2 = BCFastq.BarcodeSort(FamFq2)
    BConsFastq1, BConsFastq2 = BCFastq.pairedFastqConsolidate(BSortFq1,
                                                              BSortFq2,
                                                              stringency=0.75,
                                                              numpy=True)
    # Assuming that no reads are failed (numpy
    # consolidation does not fail reads or read pairs unless
    # there is less than 50% agreement or there are too many members
    # in a family.), just tag them, no need to check for shared pairs.
    return BConsFastq1, BConsFastq2, barcodeIndex


def pairedVCFProc(consMergeSortBAM,
                  ref="default",
                  opts="",
                  bed="default",
                  minMQ=10,
                  minBQ=20,
                  MakePileupTsv=False,
                  MakeVCF=True,
                  MakeCoverageBed=True,
                  commandStr="default"):
    if(bed == "default"):
        raise ValueError("Bed file location must be set!")
    if(ref == "default"):
        raise ValueError("Reference index location must be set!")
    # Consolidating families into single reads
    # Variant Calling Step using MPileup
    # print("Now filtering for reads with NM > 0 only if you want to.")
    Results = {}
    if(MakeCoverageBed is True):
        OutBed = PileupUtils.CalcWithinBedCoverage(consMergeSortBAM,
                                                   bed=bed,
                                                   minMQ=minMQ,
                                                   minBQ=minBQ)
    if(MakePileupTsv is True):
        PileupTSV = PileupUtils.CustomPileupToTsv(consMergeSortBAM,
                                                  bedfile=bed,
                                                  minMQ=minMQ,
                                                  minBQ=minBQ)
        pl("PileupTSV: {}".format(PileupTSV))
        Results["tsv"] = PileupTSV
    if(MakeVCF is True):
        SNP_VCF = VCFWriters.SNVCrawler(consMergeSortBAM,
                                        minMQ=minMQ,
                                        minBQ=minBQ,
                                        reference=ref,
                                        commandStr=commandStr,
                                        reference_is_path=True,
                                        bed=bed)
        CleanedVCF = BCVCF.FilterVCFFileByBed(SNP_VCF, bed)
        pl("SNP VCF: {}".format(SNP_VCF))
        Results["vcf"] = SNP_VCF
        VCFStatsFile = VCFStats(SNP_VCF)
    # AlleleFreqTSV = PileupUtils.AlleleFrequenciesByBase(consMergeSortBAM,
    #                                                     bedfile=bed)
    # This is probably useless given that I'm doing this "manually",
    # but I'm keeping this in here for good measure.
    return Results
