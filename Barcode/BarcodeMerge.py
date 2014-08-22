#!/mounts/anaconda/bin/python

from Bio import SeqIO
import argparse
import pysam
import BarcodeBamtools
import BarcodeVCFTools
import BarcodeUtils
import BarcodeHTSTools
import BarcodeFastqTools
import re
import subprocess
#Contains utilities for the completion of a variety of 
#tasks related to barcoded protocols for ultra-low
#frequency variant detection, particularly for circulating tumor DNA
#

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--fq', help="Provide your fastq file(s). Note: if '--BAM' option is set, this needs to be the Family Fastq files", nargs = "+", metavar = ('reads'),required=True)
    parser.add_argument('-r','--ref',help="Prefix for the reference index. Required.",required=True)
    parser.add_argument('--adapter', help="Adapter for samples. If not set, defaults to CAGT.",metavar=('Adapter'),default="CAGT")
    parser.add_argument('-l','--bar-len', help="Length of anticipated barcode. Defaults to 12.", default=12, type=int)
    parser.add_argument('-p','--paired-end',help="Whether the experiment is paired-end or not. Default: True",default=True)
    parser.add_argument('-k','--keepFailed', help = "To keep reads which fail filters, but leave them marked. Default: True", default=True)
    parser.add_argument('-a','--aligner', help="Provide your aligner. E.g., 'bwa', 'bowtie2', or 'snap'. Currently only 'bwa' is supported. Default: bwa", nargs='?', metavar='aligner', default='bwa')
    parser.add_argument('-o','--opts', help="Please place additional arguments for the aligner after this tag in quotation marks. E.g.: --opts '-L 0' ", nargs='?', default='')
    parser.add_argument('-b','--BAM', help="BAM file, if alignment has already run and fastq family files created.",default="default")
    parser.add_argument('-s','--sam-file',help="Name for intermediate SAM file.",default="default")
    parser.add_argument('--bed',help="full path to bed file used for variant-calling steps.",default="/yggdrasil/workspace/Barcode_Noah/WorkDir/xgen-pan-cancer-sorted.trim.bed")
    args=parser.parse_args()
    adapter=args.adapter
    
    #If single-end
    if(len(args.fq)==1):
        if(args.BAM=="default"):
            if(args.paired_end==True or args.paired_end=="True"):
                raise NameError("You provided only one fastq, but you indicated that it was paired-end. Try again!")
            if(args.keepFailed==False):
                print("For some reason, I am using an old protocol which fails to take advantage of paired-end libraries.")
                Regex1,Regex2,Hits,Misses=BarcodeFastqTools.FastqRegex(args.fq[0],adapter)
                print("Regex operation complete. Now locating adapter sequence.")
            else:
                print("Regex operation avoided for compatibility.")
                Hits = args.fq[0]
            StdFilenames,ElseFilenames=BarcodeFastqTools.AdapterLoc(Hits,adapter=adapter,keepFailed=args.keepFailed)
            print("Adapter sequences located. Hits and misses (correct location vs. incorrect location or not found) parsed out.")
            print("Now removing the adapter and the barcode.")
            tags, trimfq = BarcodeFastqTools.TrimAdapter(StdFilenames,adapter)
            print("Now generating the barcode index.")
            BarcodeIndex = BarcodeFastqTools.GenerateSingleBarcodeIndex(tags)
            FamilyFastq,TotalReads,ReadsWithFamilies = BarcodeFastqTools.GetFamilySizeSingle(trimfq,BarcodeIndex,keepFailed=args.keepFailed)
            outsam=args.sam_file
            if(args.sam_file=="default"):
                outsam = '.'.join(trimfq.split('.')[0:-1])+'.'+args.aligner+'.sam'
            outbam = '.'.join(outsam.split('.')[0:-1])+'.bam'
            print("The output SAM file with be {}, while the output BAM file will be {}".format(outsam,outbam))
            if(args.aligner=="bwa"):
                bwa_command = BarcodeHTSTools.align_bwa_se(FamilyFastq,args.ref,args.opts,outsam)
                print("Aligner command was {}".format(bwa_command))
            else:
                raise NameError("Sorry, I haven't bothered to handle other aligners yet. Whoops! Remember, I warned you about this in the help menu")
            print("Converting SAM to BAM")
            BarcodeBamtools.Sam2Bam(outsam, outbam)
        else:
            outbam=args.BAM
            FamilyFastq=args.fq[0]
        print("Now tagging reads with barcodes, family counts, and a pass/fail for the \"adapter sequence\" (not really that, but it's a working appellation).")
        taggedBAM = BarcodeBamtools.singleBarcodeTagging(FamilyFastq,outbam)
        print("Now filtering based on complexity of barcodes, the adapter presence, and a reasonably-sized family.")
        mappedPassingBarcodes,failures = BarcodeBamtools.singleFilterBam(taggedBAM,criteria="complexity,adapter,family") #Barcodes must be the same on pairs, no homopolymers of >=10, adapter must be found in the correct location
        print("Consolidating families.")
        ConsBamSingle = BarcodeBamtools.SingleConsolidate(mappedPassingBarcodes)
        print("Converting BAM to fastq.")
        consulateFastq = BarcodeBamtools.SamtoolsBam2fq(ConsBamSingle, '.'.join(ConsBamSingle.split('.')[0:-1]) + ".cons.fastq")
        print("Realigning.")
        bwa_commandCons = BarcodeHTSTools.align_bwa_se(consulateFastq,args.ref,args.opts,'.'.join( consulateFastq.split('.')[0:-1] ) )
        print("Removing reads which aren't different from the reference.")
        dissentingCons,boringCons = BarcodeBamtools.singleFilterBam(ConsBamSingle,criteria='editdistance')
        print("Now sorting reads by coordinate to prepare for MPileup.")
        CorrCons = BarcodeBamtools.CorrSort(dissentingCons)
        print("Now creating a VCF using mpileup for variant calling.")
        MPileupVCF = BarcodeVCFTools.MPileup(CorrCons, args.bed, args.ref)
        print("Initial mpileup VCF is at {}. Now removing entries which have no information.".format(MPileupVCF))
        ParsedVCF = BarcodeVCFTools.ParseVCF(MPileupVCF)
        print("Now removing those entries and parsing in the VCF Data")
        ParsedVCF.cleanRecords() #Removes entries in the VCF where there is no variant
        print("Now removing entries from the VCF and writing to a new file.")
        CleanParsedVCF = BarcodeVCFTools.CleanupPileup(MPileupVCF)
        print("This is far as the program goes at this point. Thank you for playing!")
        return
    
    #If paired-end
    elif(len(args.fq)==2):
        if(args.BAM=="default"):
            #Section 1: Completes BarcodeUtils processing of FASTQ files with barcodes
            if(args.paired_end=="False" or args.paired_end==False):
                raise NameError("You provided two fastq files, but you did not select paired-end as an option. What's up with that?")
            print("Regex operation avoided for compatibility.")
            Hits1 = args.fq[0]
            StdFilenames1,ElseFilenames1=BarcodeFastqTools.AdapterLoc(Hits1,adapter=adapter,keepFailed=args.keepFailed)
            print("Adapter sequences located. Hits and misses (correct location vs. incorrect location or not found) parsed out.")
            print("Now removing the adapter and the barcode.")
            tags1, trimfq1 = BarcodeFastqTools.TrimAdapter(StdFilenames1,adapter)
            print("Now generating the barcode index.")
            
            #For end 2
            print("Regex operation avoided for compatibility.")
            Hits2 = args.fq[1]
            StdFilenames2,ElseFilenames2=BarcodeFastqTools.AdapterLoc(Hits2,adapter=adapter,keepFailed=args.keepFailed)
            print("Adapter sequences located. Hits and misses (correct location vs. incorrect location or not found) parsed out.")
            print("Now removing the adapter and the barcode.")
            tags2, trimfq2 = BarcodeFastqTools.TrimAdapter(StdFilenames2,adapter)
            
            #Section 2: Completes Alignment
            print("Now commencing paired-end alignment with extra options: {}.".format(args.opts))
            outsam=args.sam_file
            if(args.sam_file=="default"):
                outsam = '.'.join(trimfq1.split('.')[0:-1])+'.'+args.aligner+'.sam'
            outbam = '.'.join(outsam.split('.')[0:-1])+'.bam'
            print("The output SAM file with be {}, while the output BAM file will be {}".format(outsam,outbam))
            if(args.aligner=="bwa"):
                bwa_command = BarcodeHTSTools.align_bwa(trimfq1,trimfq2,args.ref,args.opts,outsam)
                print("Aligner command was {}".format(bwa_command))
            else:
                raise ValueError("Sorry, I haven't bothered to handle other aligners yet. Whoops! Remember, I warned you about this in the help menu")
            #Section 3:Completes SAM tagging
            print("Converting SAM to BAM")
            BarcodeBamtools.Sam2Bam(outsam, outbam)
        else:
            outbam=args.BAM
            trimfq1=args.fq[0]
            trimfq2=args.fq[1]
        #TODO: Remove reads in the region of interest
        print("Now tagging reads with barcodes, family counts, and a pass/fail for the \"adapter sequence\" (not really that, but it's a working appellation).")
        taggedBAM = BarcodeBamtools.pairedBarcodeTagging(trimfq1,trimfq2,outbam)
        print("Now splitting the BAM into read 1 and read 2 files.")
        read1BAM, read2BAM = BarcodeBamtools.splitBAMByReads(taggedBAM)
        print("Now merging the barcodes from both BAM files and writing to a BAM file with a BS (Barcode Sequence) tag of the barcodes from the read and its mate.")
        concatBS = BarcodeBamtools.mergeBarcodes(read1BAM,read2BAM)
        print("BAM with merged barcodes is {}".format(concatBS))
        print("Now generating double barcode index.")
        mappedPassingBarcodes,failures = BarcodeBamtools.pairedFilterBam(concatBS,criteria="complexity,adapter,barcode") #Barcodes must be the same on pairs, no homopolymers of >=10, adapter must be found in the correct location
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
        
        #Consolidating families into single reads
        consolidatedFamilies = BarcodeBamtools.Consolidate(sortedByBarcode)
        consulate1,consulate2 = BarcodeBamtools.splitBAMByReads(consolidatedFamilies)
        consulateFastq1 = BarcodeBamtools.SamtoolsBam2fq(consulate1, '.'.join(consulate1.split('.')[0:-1]) + ".cons.fastq")
        consulateFastq2 = BarcodeBamtools.SamtoolsBam2fq(consulate2, '.'.join(consulate2.split('.')[0:-1]) + ".cons.fastq")
        consSam = '.'.join(consulate1.split('.')[0:-1]) + "cons.bwa.sam"
        consBam = '.'.join(consulate1.split('.')[0:-1]) + "cons.bwa.bam"
        bwa_command = BarcodeHTSTools.align_bwa(consulateFastq1,consulateFastq2,args.ref,args.opts,consSam)
        commandStr, consBam = BarcodeBamtools.Sam2Bam(consSam,consBam) 

        ####Variant Calling Step using MPileup
        print("Now filtering for reads with NM > 0")
        dissentingCons,boringCons = BarcodeBamtools.pairedFilterBam(consBam,criteria='editdistance')
        print("Dissenting consolidated families are in {}, while the mindless meat puppets are in {}".format(dissentingCons,boringCons))
        print("Now sorting reads by coordinate to prepare for MPileup.")
        CorrCons = BarcodeBamtools.CorrSort(dissentingCons)
        
        print("Now creating a VCF using mpileup for variant calling.")
        MPileupVCF = BarcodeVCFTools.MPileup(CorrCons, args.bed, args.ref)
        print("Initial mpileup VCF is at {}. Now removing entries which have no information.".format(MPileupVCF))
        ParsedVCF = BarcodeVCFTools.ParseVCF(MPileupVCF)
        ParsedVCF.cleanRecords() #Removes entries in the VCF where there is no variant
        CleanParsedVCF = BarcodeVCFTools.CleanupPileup(MPileupVCF)
        
        print("This is far as the program goes at this point. Thank you for playing!")
        return
    else:
        raise NameError("0k4y, sm4r7 guy - what's ^ with providing me >2 fastq files??/??")

if(__name__=="__main__"):
    main()
