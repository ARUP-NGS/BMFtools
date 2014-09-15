import logging
import re
import subprocess

import BarcodeBamtools
import BarcodeFastqTools
import BarcodeHTSTools
import BarcodeVCFTools


def pairedBamProc(consfq1, consfq2, consfqSingle="default", opts="",
                  bamPrefix="default", ref="default", aligner="default"):
    if(ref == "default"):
        raise ValueError("Reference index required!")
    if(bamPrefix == "default"):
        bamPrefix = '.'.join(consfq1.split('.')[0:-1])
    if(aligner == "default"):
        logging.info("No aligner set, defaulting to bwa.")
        aligner = "bwa"
    outsamProperPair = bamPrefix + '.sam'
    outbamProperPair = bamPrefix + '.bam'
    outsamSingle = bamPrefix + "solo.sam"
    outbamSingle = bamPrefix + "solo.bam"
    logging.info("The output SAM file: {}. Output BAM file: {}".format(
                 outsamProperPair, outbamProperPair))
    if(aligner == "bwa"):
        outsamProperPair, bwa_command = BarcodeHTSTools.align_bwa(
            consfq1, consfq2, ref, opts, outsamProperPair)
        logging.info("Aligner command was {}".format(bwa_command))
        if(consfqSingle != "default"):
            outsamSingle, bwase_command = BarcodeHTSTools.align_bwa_se(
                consfqSingle, ref, opts, outsamSingle)
            logging.info(
                "Aligner command for single-end was {}".format(bwase_command))
            logging.info("Converting single-end sam to bam")
            BarcodeBamtools.Sam2Bam(outsamSingle, outbamSingle)
            logging.info("Tagging solo BAM")
            taggedSingleBAM = BarcodeBamtools.singleBarcodeTagging(
                consfqSingle, outbamSingle)
            logging.info("Removing unmapped reads and those failing filters.")
            passTaggedSingleBAM = BarcodeBamtools.singleFilterBam(
                taggedSingleBAM,
                criteria="complexity,adapter,ismapped")
            logging.info("Now generating index for solo reads")
            singleIndex = BarcodeBamtools.GenerateBarcodeIndexBAM(
                passTaggedSingleBAM)
            logging.info("Now tagging BAM file with family size.")
            familySizeSoloBAM = BarcodeBamtools.getFamilySizeBAM(
                passTaggedSingleBAM, singleIndex)
            sortFSSBam = BarcodeBamtools.CorrSort(familySizeSoloBAM)
    else:
        raise ValueError("Sorry, only bwa is supported currently.")
    logging.info("Converting SAM to BAM")
    BarcodeBamtools.Sam2Bam(outsamProperPair, outbamProperPair)
    logging.info("Now tagging BAM with custom SAM tags.")
    taggedBAM = BarcodeBamtools.pairedBarcodeTagging(
        consfq1, consfq2, outbamProperPair)
    logging.info("Now splitting the BAM into read 1 and read 2 files.")
    read1BAM, read2BAM = BarcodeBamtools.splitBAMByReads(taggedBAM)
    logging.info("Now merging the barcodes from each pair in the BAM files.")
    concatBS = BarcodeBamtools.mergeBarcodes(read1BAM, read2BAM)
    logging.info("BAM with merged barcodes is {}".format(concatBS))
    logging.info("Now generating double barcode index.")
    mappedMerge, failures = BarcodeBamtools.pairedFilterBam(
        concatBS, criteria="complexity,adapter,barcode")
    doubleIndex = BarcodeBamtools.GenerateBarcodeIndexBAM(mappedMerge)
    p = subprocess.Popen(["wc", "-l", doubleIndex], stdout=subprocess.PIPE)
    out, err = p.communicate()
    logging.info("Number of families found: {}".format(
        re.findall(r'\d+', out)[0]))
    histochart = BarcodeBamtools.GenerateFamilyHistochart(doubleIndex)
    logging.info("Histochart of family sizes: {}".format(histochart))
    #UNCOMMENT THIS BLOCK IF YOU WANT TO START MESSING WITH RESCUE
    '''
        logging.info("Rescue step, marking the BD as their Hamming distance.")
        newRef = GenerateBarcodeIndexReference(uniqueBigFamilies)
        indexBowtie(newRef)
        mergedFastq = mergeSequencesFastq(tags1, tags2,)
        joiningSAM = CustomRefBowtiePaired(mergedFastq,newRef)
        return
        joinedFamilies = fuzzyJoining(familyMarked,joiningSAM)
        logging.info("joinedFamilies is {}".format(joinedFamilies))
    '''
    logging.info("Now determining family size for the doubled barcodes.")
    families, BCList = BarcodeBamtools.getFamilySizeBAM(
        mappedMerge, doubleIndex)
    familyP, familyF = BarcodeBamtools.pairedFilterBam(
        families, criteria="family")
    corrSorted = BarcodeBamtools.CorrSort(familyP)
    if(consfqSingle != "default"):
        mergedSinglePair = BarcodeBamtools.mergeBams(corrSorted, sortFSSBam)
        return mergedSinglePair
    return corrSorted


def pairedFastqProc(inFastq1, inFastq2, homing="default"):
    if(homing == "default"):
        homing = "CAGT"
    #For reads 1
    homingP1, homingF1 = BarcodeFastqTools.HomingSeqLoc(
        inFastq1, homing=homing)
    logging.info("Homing sequences located, reads parsed out.")
    logging.info("Now removing the homing sequence and the barcode.")
    tags1, trimfq1 = BarcodeFastqTools.TrimHoming(homingP1, homing)
    logging.info("Now generating the barcode index.")
    BarcodeIndex1 = BarcodeFastqTools.GenerateSingleBarcodeIndex(tags1)
    FamFq1, AllRds1, FamRds1 = BarcodeFastqTools.GetFamilySizeSingle(trimfq1, BarcodeIndex1)
    BarcodeSortedFastq1 = BarcodeFastqTools.BarcodeSort(FamFq1)
    #For reads 2
    homingP2, homingF2 = BarcodeFastqTools.HomingSeqLoc(
        inFastq2, homing=homing)
    logging.info("Homing sequences located, parsing reads.")
    logging.info("Now removing the homing sequence and the barcode.")
    tags2, trimfq2 = BarcodeFastqTools.TrimHoming(homingP2, homing)
    logging.info("Now generating the barcode index.")
    BarcodeIndex2 = BarcodeFastqTools.GenerateSingleBarcodeIndex(tags2)
    FamilyFastq2, TotalReads2, ReadsWithFamilies2 = BarcodeFastqTools.GetFamilySizeSingle(trimfq2, BarcodeIndex2)
    BarcodeSortedFastq2 = BarcodeFastqTools.BarcodeSort(FamilyFastq2)
    BarcodeConsFastq1, BarcodeConsFastq2 = BarcodeFastqTools.pairedFastqConsolidate(BarcodeSortedFastq1, BarcodeSortedFastq2, stringency=0.75)
    BarcodeConsFqIndex1 = BarcodeFastqTools.GenerateFullFastqBarcodeIndex(BarcodeConsFastq1)
    BarcodeConsFqIndex2 = BarcodeFastqTools.GenerateFullFastqBarcodeIndex(BarcodeConsFastq2)
    BarcodeConsPair1, BarcodeConsPair2, BarcodeSingle = BarcodeFastqTools.findProperPairs(BarcodeConsFastq1, BarcodeConsFastq2, index1=BarcodeConsFqIndex1, index2=BarcodeConsFqIndex2)
    return BarcodeConsPair1, BarcodeConsPair2, BarcodeSingle


def pairedVCFProc(consMergeSortBAM, ref="", opts="", bed=""):
    if(bed == ""):
        raise ValueError("Bed file location must be set!")
    if(ref == ""):
        raise ValueError("Reference index location must be set!")
    #Consolidating families into single reads
    ####Variant Calling Step using MPileup
    #print("Now filtering for reads with NM > 0 only if you want to.")
    logging.info("Now sorting reads by coordinate to prepare for MPileup.")
    logging.info("Now creating a VCF using mpileup for variant calling.")
    MPileupVCF = BarcodeVCFTools.MPileup(consMergeSortBAM, ref, bed=bed)
    logging.info("Initial mpileup VCF is at {}. Now removing entries which have no information.".format(MPileupVCF))
    ParsedVCF = BarcodeVCFTools.ParseVCF(MPileupVCF)
    ParsedVCF.cleanRecords()  # Removes entries in the VCF where there is no variant
    CleanParsedVCF = BarcodeVCFTools.CleanupPileup(MPileupVCF)
    return CleanParsedVCF


def singleBamProc(FamilyFastq, ref, opts, aligner="bwa", bamPrefix="default"):
    logging.info("Now tagging reads with barcodes, family counts, and a pass/fail for the homing sequence.")
    logging.info("Now filtering based on complexity of barcodes, the homing presence, and a reasonably-sized family.")
    if(bamPrefix == "default"):
        bamPrefix = FamilyFastq.split('.')[0] + '.FMS'
        outsam, outbam = bamPrefix + '.sam', bamPrefix + '.bam'
    logging.info("The output SAM file with be {}, while the output BAM file will be {}".format(outsam, outbam))
    if(aligner == "bwa"):
        outsamFile, bwa_command = BarcodeHTSTools.align_bwa_se(FamilyFastq, ref, opts, outsam)
        logging.info("Aligner command was {}".format(bwa_command))
    else:
        raise BarcodeHTSTools.IllegalArgumentError("Sorry, I haven't bothered to handle other aligners yet. Whoops! Remember, I warned you about this in the help menu")
    logging.info("Converting SAM to BAM")
    BarcodeBamtools.Sam2Bam(outsam, outbam)
    taggedBAM = BarcodeBamtools.singleBarcodeTagging(FamilyFastq, outbam)
    return taggedBAM


def singleFastqProc(inFastq, homing="default"):
    if(homing == "default"):
        homing = "CAGT"
    StdFilenames,ElseFilenames = BarcodeFastqTools.HomingSeqLoc(inFastq, homing=homing)
    logging.info("Adapter sequences located. Hits and misses (correct location vs. incorrect location or not found) parsed out.")
    logging.info("Now removing the homing and the barcode.")
    tags, trimfq = BarcodeFastqTools.TrimHoming(StdFilenames, homing)
    logging.info("Now generating the barcode index.")
    BarcodeIndex = BarcodeFastqTools.GenerateSingleBarcodeIndex(tags)
    FamilyFastq, TotalReads, ReadsWithFamilies = BarcodeFastqTools.GetFamilySizeSingle(trimfq, BarcodeIndex)
    BarcodeSortedFastq = BarcodeFastqTools.BarcodeSort(FamilyFastq)
    BarcodeConsFastq = BarcodeFastqTools.singleFastqConsolidate(BarcodeSortedFastq, stringency=0.667)
    return BarcodeConsFastq

def singleVCFProc(ConsensusBam, bed, ref):
    logging.info("Now sorting reads by coordinate to prepare for MPileup.")
    CorrCons = BarcodeBamtools.CorrSort(ConsensusBam)
    logging.info("Now creating a VCF using mpileup for variant calling.")
    MPileupVCF = BarcodeVCFTools.MPileup(CorrCons, ref, bed=bed)
    logging.info("Initial mpileup VCF is at {}. Now removing entries which have no information.".format(MPileupVCF))
    ParsedVCF = BarcodeVCFTools.ParseVCF(MPileupVCF)
    logging.info("Now removing those entries and parsing in the VCF Data")
    ParsedVCF.cleanRecords() # Removes entries in the VCF where there is no variant
    logging.info("Now removing entries from the VCF and writing to a new file.")
    CleanParsedVCF = BarcodeVCFTools.CleanupPileup(MPileupVCF)
    return CleanParsedVCF
