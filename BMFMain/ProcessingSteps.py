import re
import subprocess

from MawCluster import BCBam, VCFWriters
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
                  abrapath="default"):
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
        pl("No aligner set, defaulting to bwa.")
        aligner = "bwa"
    outsamProperPair = bamPrefix + '.sam'
    outbamProperPair = bamPrefix + '.bam'
    outbamSingle = bamPrefix + "solo.bam"
    pl("The output SAM file: {}. Output BAM file: {}".format(
        outsamProperPair, outbamProperPair))
    if(aligner == "bwa"):
        outbamProperPair = HTSUtils.align_bwa(
            consfq1, consfq2, ref, opts, outsamProperPair)
        if(consfqSingle != "default"):
            outbamSingle, bwase_command = HTSUtils.align_bwa_se(
                consfqSingle, ref, opts, outbamSingle)
            pl(
                "Aligner command for single-end was {}".format(bwase_command))
            pl("Tagging solo BAM")
            taggedSingleBAM = BCBam.singleBarcodeTagging(
                consfqSingle, outbamSingle)
            pl("Removing unmapped reads and those failing filters.")
            passTaggedSingleBAM, failTSB = BCBam.singleFilterBam(
                taggedSingleBAM,
                criteria="complexity,adapter,ismapped")
            pl("Now tagging BAM file with family size.")
            familySizeSoloBAM, famLst = BCBam.getFamilySizeBAM(
                passTaggedSingleBAM, barIndex)
            sortFSSBam = BCBam.CoorSort(familySizeSoloBAM)
    else:
        raise ValueError("Sorry, only bwa is supported currently.")
    pl("Now tagging BAM with custom SAM tags.")
    taggedBAM = BCBam.pairedBarcodeTagging(
        consfq1, consfq2, outbamProperPair)
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
    mappedPass, failures = BCBam.pairedFilterBam(
        namesortedRealignedFull, criteria="adapter,ismapped")
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
    pl("Now determining family size for the doubled barcodes.")
    families, BCList = BCBam.getFamilySizeBAM(
        mappedPass, barIndex)
    # This step not needed for shades protocol, as fastq
    # families have already been filtered for size.
    # familyP, familyF = BCBam.pairedFilterBam(
    #    families, criteria="family")
    SVBam, MarkedFamilies = SVRP(families,
                                 bedfile=bed,
                                 tempBAMPrefix=families[0:-4],
                                 summary=(families[0:-4] +
                                          '.SV.txt'))
    pl(("{} is the bam with all reads considered relevant ".format(SVBam) +
        "to translocations."))
    # SVOutputFile = BCBam.CallTranslocations(SVBam, bedfile=bed)
    coorSorted = BCBam.CoorSort(MarkedFamilies)
    CoverageBed = PileupUtils.BamToCoverageBed(coorSorted, mincov=mincov)
    pl("Coverage bed: {}".format(CoverageBed))
    if(consfqSingle != "default"):
        mergedSinglePair = BCBam.mergeBams(coorSorted, sortFSSBam)
        return mergedSinglePair
    return coorSorted


def pairedFastqShades(inFastq1, inFastq2, indexFastq, stringency=0.75):
    bcFastq1, bcFastq2 = BCFastq.FastqPairedShading(inFastq1,
                                                    inFastq2,
                                                    indexFastq,
                                                    gzip=False)
    pl("Beginning pairedFastqShades for {}, {}".format(inFastq1, inFastq2))
    barcodeIndex = BCFastq.GenerateShadesIndex(indexFastq)
    (FamFqs, SingleFqs, numReads,
     numReadsWFam) = BCFastq.GetFamilySizePaired(bcFastq1,
                                                 bcFastq2, barcodeIndex)
    FamFq1 = FamFqs[0]
    FamFq2 = FamFqs[1]
    SingleFq1 = SingleFqs[0]
    SingleFq2 = SingleFqs[1]
    """
    TODO: Write this step
    NumberRescued = BCFastq.ShadesRescuePaired(SingleFq1, SingleFq2,
                                               appendFq1=FamFq1,
                                               appendFq2=FamFq2,
                                               index=barcodeIndex)
    """
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
                  reference="default",
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
                                        reference=reference,
                                        commandStr=commandStr,
                                        reference_is_path=True,
                                        bedfile=bed)
        pl("SNP VCF: {}".format(SNP_VCF))
        Results["vcf"] = SNP_VCF
        VCFStatsFile = VCFStats(SNP_VCF)
    # AlleleFreqTSV = PileupUtils.AlleleFrequenciesByBase(consMergeSortBAM,
    #                                                     bedfile=bed)
    # This is probably useless given that I'm doing this "manually",
    # but I'm keeping this in here for good measure.
    return Results
