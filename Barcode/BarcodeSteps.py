#/mounts/anaconda/bin/python

import subprocess
import re

import BarcodeBamtools
import BarcodeFastqTools
import BarcodeHTSTools
import BarcodeVCFTools

def pairedBamProc(trimfq1, trimfq2,outbam):
    print("Now tagging reads with barcodes, family counts, and a pass/fail for the homing sequence.")
    taggedBAM = BarcodeBamtools.pairedBarcodeTagging(trimfq1,trimfq2,outbam)
    print("Now splitting the BAM into read 1 and read 2 files.")
    read1BAM, read2BAM = BarcodeBamtools.splitBAMByReads(taggedBAM)
    print("Now merging the barcodes from both BAM files and writing to a BAM file with a BS (Barcode Sequence) tag of the barcodes from the read and its mate.")
    concatBS = BarcodeBamtools.mergeBarcodes(read1BAM,read2BAM)
    print("BAM with merged barcodes is {}".format(concatBS))
    print("Now generating double barcode index.")
    mappedPassingBarcodes,failures = BarcodeBamtools.pairedFilterBam(concatBS,criteria="complexity,adapter,barcode") #Barcodes must be the same on pairs, no homopolymers of >=10, homing must be found in the correct location
    doubleIndex = BarcodeBamtools.GenerateBarcodeIndexBAM(mappedPassingBarcodes)
    p = subprocess.Popen(["wc","-l",doubleIndex],stdout=subprocess.PIPE)
    out, err = p.communicate()
    print("Number of families found: {}".format(re.findall(r'\d+',out)[0]))
    histochart = BarcodeBamtools.GenerateFamilyHistochart(doubleIndex)
    print("Histochart of family sizes available at: {}".format(histochart))
    #UNCOMMENT THIS BLOCK IF YOU WANT TO START MESSING WITH RESCUE
    '''
        print("Now bringing those within 2 base pair differences into the main group, marking the BD as their Hamming distance.")
        newRef = BarcodeBamtools.BarcodeFastqTools.GenerateBarcodeIndexReference(uniqueBigFamilies)
        indexBowtie(newRef)
        mergedFastq = mergeSequencesFastq(tags1, tags2,)
        joiningSAM = CustomRefBowtiePaired(mergedFastq,newRef)
        return
        joinedFamilies = fuzzyJoining(familyMarked,joiningSAM)
        print("joinedFamilies is {}".format(joinedFamilies))
    '''
    print("Now determining family size for the doubled barcodes.")
    families,BarcodeFamilyList = BarcodeBamtools.getFamilySizeBAM(mappedPassingBarcodes, doubleIndex)
    familyPass,familyFail = BarcodeBamtools.pairedFilterBam(families,criteria="family")
    sortedByBarcode = BarcodeBamtools.BarcodeSort(familyPass)
    return sortedByBarcode

def pairedFastqProcAlign(inFastq1,inFastq2,homing,opts="",bamPrefix="default",ref="default", aligner="default"):
    if(ref=="default"):
        raise ValueError("Reference index required!")
    if(homing == "default"):
        homing = "CAGT"
    if(bamPrefix == "default"):
        bamPrefix = '.'.join(inFastq1.split('.')[0:-1])
    if(aligner=="default"):
        aligner="bwa"
    outsam = bamPrefix + '.sam'
    outbam = bamPrefix + '.bam'
    #For reads 1
    StdFilenames1,ElseFilenames1=BarcodeFastqTools.HomingSeqLoc(inFastq1,homing=homing)
    print("Adapter sequences located. Hits and misses (correct location vs. incorrect location or not found) parsed out.")
    print("Now removing the homing sequence and the barcode.")
    tags1, trimfq1 = BarcodeFastqTools.TrimHoming(StdFilenames1,homing)
    print("Now generating the barcode index.")
    
    #For reads 2
    StdFilenames2,ElseFilenames2=BarcodeFastqTools.HomingSeqLoc(inFastq2,homing=homing)
    print("Adapter sequences located. Hits and misses (correct location vs. incorrect location or not found) parsed out.")
    print("Now removing the homing sequence and the barcode.")
    tags2, trimfq2 = BarcodeFastqTools.TrimHoming(StdFilenames2,homing)
    print("Now generating the barcode index.")
    
    print("The output SAM file with be {}, while the output BAM file will be {}".format(outsam,outbam))
    if(aligner=="bwa"):
        outsam, bwa_command = BarcodeHTSTools.align_bwa(trimfq1,trimfq2,ref,opts,outsam)
        print("Aligner command was {}".format(bwa_command))
    else:
        raise ValueError("Sorry, I haven't bothered to handle other aligners yet. Whoops! Remember, I warned you about this in the help menu")
    print("Converting SAM to BAM")
    BarcodeBamtools.Sam2Bam(outsam, outbam)
    return outbam, trimfq1, trimfq2

def pairedVCFProc(sortedByBarcode,ref="",opts="",bed=""):
    if(bed==""):
        raise ValueError("Bed file location must be set!")
    if(ref==""):
        raise ValueError("Reference index location must be set!")
    #Consolidating families into single reads
    consolidatedFamilies = BarcodeBamtools.Consolidate(sortedByBarcode)
    consulate1,consulate2 = BarcodeBamtools.splitBAMByReads(consolidatedFamilies)
    consulateFastq1 = BarcodeBamtools.SamtoolsBam2fq(consulate1, '.'.join(consulate1.split('.')[0:-1]) + ".cons.fastq")
    consulateFastq2 = BarcodeBamtools.SamtoolsBam2fq(consulate2, '.'.join(consulate2.split('.')[0:-1]) + ".cons.fastq")
    consSam = '.'.join(consulate1.split('.')[0:-1]) + "cons.bwa.sam"
    consBam = '.'.join(consulate1.split('.')[0:-1]) + "cons.bwa.bam"
    consSam, bwa_command = BarcodeHTSTools.align_bwa(consulateFastq1,consulateFastq2,ref,opts,consSam)
    print("Now attempting to convert {} to BAM, with the finished file at {}".format(consSam,consBam))
    commandStr, consBam = BarcodeBamtools.Sam2Bam(consSam,consBam) 
    ####Variant Calling Step using MPileup
    #print("Now filtering for reads with NM > 0")
    #dissentingCons,boringCons = BarcodeBamtools.pairedFilterBam(consBam,criteria='editdistance')
    #print("Dissenting consolidated families are in {}, while the mindless meat puppets are in {}".format(dissentingCons,boringCons))
    print("Now sorting reads by coordinate to prepare for MPileup.")
    CorrCons = BarcodeBamtools.CorrSort(consBam)
    
    print("Now creating a VCF using mpileup for variant calling.")
    MPileupVCF = BarcodeVCFTools.MPileup(CorrCons, bed, ref)
    print("Initial mpileup VCF is at {}. Now removing entries which have no information.".format(MPileupVCF))
    ParsedVCF = BarcodeVCFTools.ParseVCF(MPileupVCF)
    ParsedVCF.cleanRecords() #Removes entries in the VCF where there is no variant
    CleanParsedVCF = BarcodeVCFTools.CleanupPileup(MPileupVCF)
    return CleanParsedVCF

def singleBamProc(FamilyFastq,outbam,ref,opts):
    print("Now tagging reads with barcodes, family counts, and a pass/fail for the homing sequence.")
    print("Now filtering based on complexity of barcodes, the homing presence, and a reasonably-sized family.")
    taggedBAM = BarcodeBamtools.singleBarcodeTagging(FamilyFastq,outbam)
    mappedPassingBarcodes,failures = BarcodeBamtools.singleFilterBam(taggedBAM,criteria="complexity,adapter,family") #Barcodes must be the same on pairs, no homopolymers of >=10, homing sequence must be found in the correct location
    print("Sorting by barcode to prepare for consolidation.")
    sortPassingBarcodes = BarcodeBamtools.BarcodeSort(mappedPassingBarcodes, paired=False)
    print("Consolidating families.")
    ConsBamSingle = BarcodeBamtools.SingleConsolidate(sortPassingBarcodes)
    print("Converting BAM to fastq.")
    consulateFastq = BarcodeBamtools.SamtoolsBam2fq(ConsBamSingle, '.'.join(ConsBamSingle.split('.')[0:-1]) + ".cons.fastq")
    print("Consulate Fastq is at {}".format(consulateFastq))
    consensuses,bwa_commandCons = BarcodeHTSTools.align_bwa_se(consulateFastq,ref,opts,consulateFastq.split('.')[0] + ".cons.realigned.sam" )
    print("Realigned, now in {}".format(consensuses))
    print("Now sorting reads by coordinate to prepare for MPileup.")
    ConsBamStr, ConsBamRealigned = BarcodeBamtools.Sam2Bam(consensuses, '.'.join(consensuses.split('.')[0:-1])+'.bam')
    CorrCons = BarcodeBamtools.CorrSort(ConsBamRealigned)
    return CorrCons

def singleFastqProcAlign(inFastq,ref,aligner,opts,homing="default",bamPrefix="default"):
    if(homing == "default"):
        homing = "CAGT"
    if(bamPrefix == "default"):
        bamPrefix = '.'.join(inFastq.split('.')[0:-1])
    StdFilenames,ElseFilenames=BarcodeFastqTools.HomingSeqLoc(inFastq,homing=homing)
    print("Adapter sequences located. Hits and misses (correct location vs. incorrect location or not found) parsed out.")
    print("Now removing the homing and the barcode.")
    tags, trimfq = BarcodeFastqTools.TrimHoming(StdFilenames,homing)
    print("Now generating the barcode index.")
    BarcodeIndex = BarcodeFastqTools.GenerateSingleBarcodeIndex(tags)
    FamilyFastq,TotalReads,ReadsWithFamilies = BarcodeFastqTools.GetFamilySizeSingle(trimfq,BarcodeIndex)
    outsam, outbam = bamPrefix + '.sam', bamPrefix + '.bam'
    print("The output SAM file with be {}, while the output BAM file will be {}".format(outsam,outbam))
    if(aligner=="bwa"):
        outsamFile, bwa_command = BarcodeHTSTools.align_bwa_se(FamilyFastq,ref,opts,outsam)
        print("Aligner command was {}".format(bwa_command))
    else:
        raise BarcodeHTSTools.IllegalArgumentError("Sorry, I haven't bothered to handle other aligners yet. Whoops! Remember, I warned you about this in the help menu")
    print("Converting SAM to BAM")
    BarcodeBamtools.Sam2Bam(outsam, outbam)
    return outbam, FamilyFastq

def singleVCFProc(ConsensusBam,bed,ref):
    print("Now creating a VCF using mpileup for variant calling.")
    MPileupVCF = BarcodeVCFTools.MPileup(ConsensusBam, bed, ref)
    print("Initial mpileup VCF is at {}. Now removing entries which have no information.".format(MPileupVCF))
    ParsedVCF = BarcodeVCFTools.ParseVCF(MPileupVCF)
    print("Now removing those entries and parsing in the VCF Data")
    ParsedVCF.cleanRecords() #Removes entries in the VCF where there is no variant
    print("Now removing entries from the VCF and writing to a new file.")
    CleanParsedVCF = BarcodeVCFTools.CleanupPileup(MPileupVCF)
    return CleanParsedVCF
