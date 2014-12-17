import re
import subprocess

import BCBam
import BCFastq
import HTSUtils
import BCVCF
from HTSUtils import printlog as pl


def pairedBamProc(consfq1, consfq2, consfqSingle="default", opts="",
                  bamPrefix="default", ref="default", aligner="default",
                  barIndex="default",
                  bed="/yggdrasil/workspace/Barcode_Noah/cfDNA_targets.bed"):
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
    outsamSingle = bamPrefix + "solo.sam"
    outbamSingle = bamPrefix + "solo.bam"
    pl("The output SAM file: {}. Output BAM file: {}".format(
        outsamProperPair, outbamProperPair))
    if(aligner == "bwa"):
        outsamProperPair = HTSUtils.align_bwa(
            consfq1, consfq2, ref, opts, outsamProperPair)
        if(consfqSingle != "default"):
            outsamSingle, bwase_command = HTSUtils.align_bwa_se(
                consfqSingle, ref, opts, outsamSingle)
            pl(
                "Aligner command for single-end was {}".format(bwase_command))
            pl("Converting single-end sam to bam")
            BCBam.Sam2Bam(outsamSingle, outbamSingle)
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
    pl("Converting SAM to BAM")
    BCBam.Sam2Bam(outsamProperPair, outbamProperPair)
    pl("Now tagging BAM with custom SAM tags.")
    taggedBAM = BCBam.pairedBarcodeTagging(
        consfq1, consfq2, outbamProperPair)
    pl("Now splitting the BAM into read 1 and read 2 files.")
    pl("Now generating double barcode index.")
    realignedFull = BCBam.AbraCadabra(taggedBAM, ref=ref, bed=bed)
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
    '''
        pl("Rescue step, marking the BD as their Hamming distance.")
        newRef = GenerateBarcodeIndexReference(uniqueBigFamilies)
        indexBowtie(newRef)
        mergedFastq = mergeSequencesFastq(tags1, tags2,)
        joiningSAM = CustomRefBowtiePaired(mergedFastq,newRef)
        return
        joinedFamilies = fuzzyJoining(familyMarked,joiningSAM)
        pl("joinedFamilies is {}".format(joinedFamilies))
    '''
    pl("Now determining family size for the doubled barcodes.")
    families, BCList = BCBam.getFamilySizeBAM(
        mappedPass, barIndex)
    familyP, familyF = BCBam.pairedFilterBam(
        families, criteria="family")
    coorSorted = BCBam.CoorSort(familyP)
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
    '''
    TODO: Write this step
    NumberRescued = BCFastq.ShadesRescuePaired(SingleFq1, SingleFq2,
                                               appendFq1=FamFq1,
                                               appendFq2=FamFq2,
                                               index=barcodeIndex)
    '''
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


def pairedFastqProc(inFastq1, inFastq2, homing="default",
                    stringency="default", useNumpy=True):
    if(stringency == "default"):
        stringency = 0.75
    if(homing == "default"):
        homing = "CAGT"
    # For reads 1
    homingP1, homingF1 = BCFastq.HomingSeqLoc(
        inFastq1, homing=homing)
    pl("Homing sequences located, reads parsed out.")
    pl("Now removing the homing sequence and the barcode.")
    tags1, trimfq1 = BCFastq.TrimHoming(homingP1, homing)
    # For reads 2
    homingP2, homingF2 = BCFastq.HomingSeqLoc(
        inFastq2, homing=homing)
    pl("Homing sequences located, parsing reads.")
    pl("Now removing the homing sequence and the barcode.")
    tags2, trimfq2 = BCFastq.TrimHoming(homingP2, homing)
    mergeTags1, mergeTags2 = BCFastq.mergeBarcodes(trimfq1, trimfq2)
    pl("Now generating the barcode index.")
    BarcodeIndex = BCFastq.PairFastqBarcodeIndex(mergeTags1, mergeTags2)
    FamFq1, AllRds1, FamRds1 = BCFastq.GetFamilySizeSingle(mergeTags1,
                                                           BarcodeIndex)
    FamFq2, AllRds2, FamRds2 = BCFastq.GetFamilySizeSingle(mergeTags2,
                                                           BarcodeIndex)
    BSortFq1 = BCFastq.BarcodeSort(FamFq1)
    BSortFq2 = BCFastq.BarcodeSort(FamFq2)
    BConsFastq1, BConsFastq2 = BCFastq.pairedFastqConsolidate(
        BSortFq1, BSortFq2, stringency=stringency, numpy=useNumpy)
    BConsFqIndex1 = BCFastq.GenerateSingleBarcodeIndex(BConsFastq1)
    BConsFqIndex2 = BCFastq.GenerateSingleBarcodeIndex(BConsFastq2)
    sharedBC = BCFastq.getSharedBC(BConsFqIndex1, BConsFqIndex2)
    BConsPair1, BConsPair2, BarcodeSingle = BCFastq.getProperPairs(
        BConsFastq1, BConsFastq2, shared=sharedBC)
    BConsBSort1 = BCFastq.BarcodeSort(BConsPair1)
    BConsBSort2 = BCFastq.BarcodeSort(BConsPair2)
    RenameFq1, RenameFq2 = BCFastq.renameReads(BConsBSort1, BConsBSort2)
    printStr = "Now returning BConsPair1 ({}), BConsPair2, ({})".format(
        RenameFq1, RenameFq2)
    printStr += ", and BarcodeSingle ({})".format(BarcodeSingle)
    pl(printStr)
    return RenameFq1, RenameFq2, BarcodeSingle, BarcodeIndex


def pairedVCFProc(consMergeSortBAM, ref="", opts="", bed=""):
    if(bed == ""):
        raise ValueError("Bed file location must be set!")
    if(ref == ""):
        raise ValueError("Reference index location must be set!")
    # Consolidating families into single reads
    # Variant Calling Step using MPileup
    # print("Now filtering for reads with NM > 0 only if you want to.")
    pl("Now sorting reads by coordinate to prepare for MPileup.")
    pl("Now creating a VCF using mpileup for variant calling.")
    MPileupVCF = BCVCF.MPileup(consMergeSortBAM, ref, bed=bed)
    pl("Initial mpileup: {}. Filtering.".format(MPileupVCF))
    ParsedVCF = BCVCF.ParseVCF(MPileupVCF)
    return ParsedVCF


def singleBamProc(FamilyFastq, ref, opts, aligner="bwa", bamPrefix="default"):
    pl("Now tagging reads.")
    pl("Now filtering reads")
    if(bamPrefix == "default"):
        bamPrefix = FamilyFastq.split('.')[0] + '.FMS'
        outsam, outbam = bamPrefix + '.sam', bamPrefix + '.bam'
    pl("Output Sam: {}. Output Bam: {}".format(outsam, outbam))
    if(aligner == "bwa"):
        outsamFile, bwa_command = HTSUtils.align_bwa_se(
            FamilyFastq, ref, opts, outsam)
        pl("Aligner command was {}".format(bwa_command))
    else:
        raise ValueError("Sorry, I don't handle that aligner.")
    pl("Converting SAM to BAM")
    BCBam.Sam2Bam(outsam, outbam)
    taggedBAM = BCBam.singleBarcodeTagging(FamilyFastq, outbam)
    return taggedBAM


def singleFastqProc(inFastq, homing="default"):
    if(homing == "default"):
        homing = "CAGT"
    StdFilenames, ElseFilenames = BCFastq.HomingSeqLoc(inFastq, homing=homing)
    pl("Homing seq located, parsing these out.")
    pl("Now removing the homing and the barcode.")
    tags, trimfq = BCFastq.TrimHoming(StdFilenames, homing)
    pl("Now generating the barcode index.")
    BarcodeIndex = BCFastq.GenerateSingleBarcodeIndex(tags)
    FamilyFastq, TotalReads, FamReads = BCFastq.GetFamilySizeSingle(
        trimfq, BarcodeIndex)
    BSortFq = BCFastq.BarcodeSort(FamilyFastq)
    BConsFastq = BCFastq.singleFastqConsolidate(BSortFq, stringency=0.667)
    return BConsFastq


def singleVCFProc(ConsensusBam, bed, ref):
    pl("Now sorting reads by coordinate to prepare for MPileup.")
    CorrCons = BCBam.CoorSort(ConsensusBam)
    pl("Now creating a VCF using mpileup for variant calling.")
    MPileupVCF = BCVCF.MPileup(CorrCons, ref, bed=bed)
    pl("Initial mpileup, with ALT entries of X removed: {}".format(MPileupVCF))
    ParsedVCF = BCVCF.ParseVCF(MPileupVCF)
    return ParsedVCF
