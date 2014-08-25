#/mounts/anaconda/bin/python

import BarcodeBamtools
import BarcodeFastqTools
import BarcodeHTSTools
import BarcodeVCFTools

def pairedBamProc():
    return

def pairedFastqProc():
    return

def pairedVCFProc():
    return

def singleBamProc(FamilyFastq,outbam,ref,opts):
    print("Now tagging reads with barcodes, family counts, and a pass/fail for the \"adapter sequence\" (not really that, but it's a working appellation).")
    print("Now filtering based on complexity of barcodes, the adapter presence, and a reasonably-sized family.")
    taggedBAM = BarcodeBamtools.singleBarcodeTagging(FamilyFastq,outbam)
    mappedPassingBarcodes,failures = BarcodeBamtools.singleFilterBam(taggedBAM,criteria="complexity,adapter,family") #Barcodes must be the same on pairs, no homopolymers of >=10, adapter must be found in the correct location
    print("Consolidating families.")
    ConsBamSingle = BarcodeBamtools.SingleConsolidate(mappedPassingBarcodes)
    print("Converting BAM to fastq.")
    consulateFastq = BarcodeBamtools.SamtoolsBam2fq(ConsBamSingle, '.'.join(ConsBamSingle.split('.')[0:-1]) + ".cons.fastq")
    print("Realigning.")
    consensuses,bwa_commandCons = BarcodeHTSTools.align_bwa_se(consulateFastq,ref,opts,'.'.join( consulateFastq.split('.')[0:-1] ) )
    print("Now sorting reads by coordinate to prepare for MPileup.")
    CorrCons = BarcodeBamtools.CorrSort(consensuses)
    return CorrCons

def singleFastqProcAlign(inFastq,ref,homing="default",bamPrefix="default",aligner,opts):
    if(homing == "default"):
        homing = "CAGT"
    if(bamPrefix == "default"):
        bamPrefix = '.'.join(inFastq.split('.')[0:-1])
    StdFilenames,ElseFilenames=BarcodeFastqTools.AdapterLoc(inFastq,homing=homing,keepFailed=True)
    print("Adapter sequences located. Hits and misses (correct location vs. incorrect location or not found) parsed out.")
    print("Now removing the homing and the barcode.")
    tags, trimfq = BarcodeFastqTools.TrimAdapter(StdFilenames,homing)
    print("Now generating the barcode index.")
    BarcodeIndex = BarcodeFastqTools.GenerateSingleBarcodeIndex(tags)
    FamilyFastq,TotalReads,ReadsWithFamilies = BarcodeFastqTools.GetFamilySizeSingle(trimfq,BarcodeIndex,keepFailed=True)
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
    print("This is far as the program goes at this point. Thank you for playing!")
    return CleanParsedVCF