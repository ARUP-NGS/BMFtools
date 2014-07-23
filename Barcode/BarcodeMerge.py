#!/mounts/anaconda/bin/python

from Bio import SeqIO
import argparse

#Contains utilities for the completion of a variety of 
#
#

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--fq', help="Provide your fastq file(s)", nargs = "+", metavar = ('reads'),required=True)
    parser.add_argument('-r','--ref',help="Prefix for the reference index. Required.",required=True)
    parser.add_argument('-b','--adapter', help="Adapter for samples. If not set, defaults to CAGT.",metavar=('Adapter'),default="CAGT")
    parser.add_argument('-l','--bar-len', help="Length of anticipated barcode. Defaults to 12.", default=12, type=int)
    parser.add_argument('-p','--paired-end',help="Whether the experiment is paired-end or not. Default: True",default=True)
    parser.add_argument('-k','--keepFailed', help = "To keep reads which fail filters, but leave them marked. Default: True", default=True)
    parser.add_argument('-a','--aligner', help="Provide your aligner. E.g., 'bwa', 'bowtie2', or 'snap'. Currently only 'bwa' is supported. Default: bwa", nargs='?', metavar='aligner', default='bwa')
    parser.add_argument('-o','--opts', help="Please place additional arguments for the aligner after this tag in quotation marks. E.g.: --opts '-L 0' ", nargs='?', default='')
    parser.add_argument('-b','--BAM', help="BAM file, if alignment has already run.")
    parser.add_argument('-s','--sam-file',help="Name for intermediate SAM file.",default="default")
    args=parser.parse_args()
    adapter=args.adapter
    
    #If single-end
    if(len(args.fq)==1):
        if(args.paired_end==True or args.paired_end=="True"):
            raise NameError("You provided only one fastq, but you indicated that it was paired-end. Try again!")
        if(args.keepFailed==False):
            print("For some reason, I am using an old protocol which fails to take advantage of paired-end libraries.")
            Regex1,Regex2,Hits,Misses=FastqRegex(args.fq[0],adapter)
            print("Regex operation complete. Now locating adapter sequence.")
        else:
            print("Regex operation avoided for compatibility.")
            Hits = args.fq[0]
        StdFilenames,ElseFilenames=AdapterLoc(Hits,adapter=adapter,keepFailed=args.keepFailed)
        print("Adapter sequences located. Hits and misses (correct location vs. incorrect location or not found) parsed out.")
        print("Now removing the adapter and the barcode.")
        tags, trimfq = TrimAdapter(StdFilenames,adapter)
        BarcodeIndex = GenerateBarcodeIndex(tags)
        FamilyFastq,TotalReads,ReadsWithFamilies = GetFamilySize(trimfq,BarcodeIndex,keepFailed=args.keepFailed)
        raise NameError("Okay, so I am being pretty lazy and didn't finish the program for single-end just yet. This fatal error is not your fault, but mine. Mea culpa.")
    #If paired-end
    elif(len(args.fq)==2):
        
        #Section 1: Completes BarcodeUtils processing of FASTQ files with barcodes
        if(args.paired_end=="False" or args.paired_end==False):
            raise NameError("You provided two fastq files, but you did not select paired-end as an option. What's up with that?")
        if(args.keepFailed==False):
            print("For some reason, I am using an old protocol which fails to take advantage of paired-end libraries.")
            Regex11,Regex21,Hits1,Misses1=FastqRegex(args.fq[0],adapter)
            print("Regex operation complete. Now locating adapter sequence.")
        else:
            print("Regex operation avoided for compatibility.")
            Hits1 = args.fq[0]
        StdFilenames1,ElseFilenames1=AdapterLoc(Hits1,adapter=adapter,keepFailed=args.keepFailed)
        print("Adapter sequences located. Hits and misses (correct location vs. incorrect location or not found) parsed out.")
        print("Now removing the adapter and the barcode.")
        tags1, trimfq1 = TrimAdapter(StdFilenames1,adapter)
        BarcodeIndex1 = GenerateBarcodeIndex(tags1)
        FamilyFastq1,TotalReads1,ReadsWithFamilies1 = GetFamilySize(trimfq1,BarcodeIndex1,keepFailed=args.keepFailed)
        print("Total number of reads in read file 1 is {}, whereas the number of reads with families is {} ".format(TotalReads1,ReadsWithFamilies1))
        if(args.keepFailed==False):
            print("For some reason, I am using an old protocol which fails to take advantage of paired-end libraries.")
            Regex12,Regex22,Hits2,Misses2=FastqRegex(args.fq[1],adapter)
            print("Regex operation complete. Now locating adapter sequence.")
        else:
            print("Regex operation avoided for compatibility.")
            Hits1 = args.fq[1]
        StdFilenames2,ElseFilenames2=AdapterLoc(Hits2,adapter=adapter,keepFailed=args.keepFailed)
        print("Adapter sequences located. Hits and misses (correct location vs. incorrect location or not found) parsed out.")
        print("Now removing the adapter and the barcode.")
        tags2, trimfq2 = TrimAdapter(StdFilenames2,adapter)
        BarcodeIndex2 = GenerateBarcodeIndex(tags2)
        FamilyFastq2,TotalReads2,ReadsWithFamilies2 = GetFamilySize(trimfq1,BarcodeIndex1,keepFailed=args.keepFailed)
        print("Total number of reads in read file 2 is {}, whereas the number of reads with families is {} ".format(TotalReads2,ReadsWithFamilies2))
        
        #Section 2: Completes Alignment
        outsam=args.sam
        if(args.sam=="default"):
            outsam = '.'.join(trimfq.split('.')[0:-1])+'.'+args.aligner+'.sam'
        outbam = '.'.join(outsam.split('.')[0:-1])+'.bam'
        if(args.aligner=="bwa"):
            bwa_command = align_bwa(FamilyFastq1,FamilyFastq2,args.ref,args.opts,outsam)
            print("Aligner command was {}".format(bwa_command))
        else:
            raise NameError("Sorry, I haven't bothered to handle other aligners yet. Whoops! Remember, I warned you about this in the help menu")
        #Section 3:Completes SAM tagging
        Sam2Bam(outsam, outbam)
        taggedBAM = pairedBarcodeBamtools(FamilyFastq1,FamilyFastq2,outbam)
        read1BAM, read2BAM = splitBAMByReads(taggedBAM)
        return
    else:
        raise NameError("0k4y, sm4rt guy - what's ^ with providing me >2 fastq files?")
    return

def splitBAMByReads(BAM,read1BAM="default",read2BAM="default"):
    presplitBAM = Samfile(BAM,"rb")
    if(read1BAM=="default"):
        read1BAM = '.'.join(BAM.split('.')[0:-1])+'.R1.bam'
    if(read2BAM=="default"):
        read2BAM = '.'.join(BAM.split('.')[0:-1])+'.R2.bam'
    out1 = Samfile(read1BAM,"wb",template=presplitBAM)
    out2 = Samfile(read2BAM,"wb",template=presplitBAM)
    for entry in presplitBAM:
        if(entry.is_read1):
            out1.write(entry)
        if(entry.is_read2):
            out2.write(entry)
    out1.close()
    out2.close()
    return 

def pairedBarcodeBamtools(fq1,fq2,bam,outputBAM="default"):
    if(outputBAM=="default"):
        outputBAM = '.'.join(bam.split('.')[0:-1]) + "tagged.bam"
    read1 = SeqIO.parse(fq1, "fastq")
    read2 = SeqIO.parse(fq2, "fastq")
    #inBAM = removeSecondary(args.bam_file) #Artefactual code
    postFilterBAM = Samfile(bam,"rb")
    output=outputBAM
    if(output=="default"):
        output='.'.join(bam.split('.')[0:-1])+'.tagged.bam'
    outBAM = Samfile(output,"wb",template=postFilterBAM)
    for entry in postFilterBAM:
        if(entry.is_secondary):
            continue
        if(entry.is_read1):
            #print("Read is read 1. Now getting variables from read 1's files for the tags.")
            tempRead = read1.next()
        elif(entry.is_read2):
            #print("Read is read 2. Now getting variables from read 2's files for the tags.")
            tempRead = read2.next()
        descArray=tempRead.description.split("###")
        entry.tags = entry.tags + [("BS",descArray[-2].strip())]
        entry.tags = entry.tags + [("FM",descArray[-1].strip())]
        if(descArray[-3].strip() == "AdapterPass"):
            entry.tags = entry.tags + [("AL",1)]
        else:
            entry.tags = entry.tags + [("AL",0)]
        outBAM.write(entry)
    outBAM.close()
    postFilterBAM.close()
    return outputBAM

def removeSecondary(inBAM,outBAM="default"):
    print("Attempting to remove secondary")
    from subprocess import call
    if(outBAM=="default"):
        outBAM = '.'.join(inBAM.split('.')[0:-1])+'.2ndrm.bam'
    command = "samtools view -hb -F 0x0100 {} -o {}".format(inBAM,outBAM)
    call(command,shell=True)
    print(command + "is shell call")
    return outBAM

def GetFamilySize(trimfq,BarcodeIndex,outfq="default",singlefq="default",keepFailed=True):
    infq = SeqIO.parse(trimfq, "fastq")
    if(outfq=="default"):
        outfq = '.'.join(trimfq.split('.')[0:-1])+".fam.fastq"
    outfqBuffers = open(outfq,"w",0)
    index = open(BarcodeIndex,"r")
    if(singlefq=="default"):
        singlefq = '.'.join(trimfq.split('.')[0:-1])+".lonely.hearts.club.band.fastq"
    singlefqBuffer = open(singlefq,"w",0)
    TotalReads, ReadsWithFamilies = 0,0
    dictEntries = [line.split() for line in index]
    BarDict = {}
    for entry in dictEntries:
        BarDict[entry[1]]=entry[0]
    for read in infq:
        index.seek(0)
        TotalReads+=1
        #print("description is {}".format(read.description)) #FOR DEBUGGING PURPOSES, UNCOMMENT
        #print("-1 index of description is {}".format(read.description.split("###")[-1].strip()))
        #print("-2 index of description is {}".format(read.description.split("###")[-2].strip()))
        readTag = read.description.split("###")[-1].strip()
        newRead = read
        #print("readTag is _" + readTag + "_")
        famSize = BarDict[readTag]
        newRead.description = read.description+" ###"+famSize
        #print("famSize = _{}_".format(str(famSize)))
        #print("The value of this comparison to 1 is {}".format(str(famSize=="1")))
        if(famSize == "1"):
            SeqIO.write(newRead,singlefqBuffer,"fastq")
            if(keepFailed==True):
                SeqIO.write(newRead,outfqBuffers,"fastq")
        else:
            SeqIO.write(newRead,outfqBuffers,"fastq")
    return outfq,TotalReads, ReadsWithFamilies

def GenerateBarcodeIndex(tags_file,index_file="default"):
    from subprocess import call
    if(index_file=="default"):
        index_file = '.'.join(tags_file.split('.')[0:-1]) + ".barIdx"
    call("cat {} | sed 's:###: ###:g' | grep -v \"AdapterFail\" | paste - - - - | awk 'BEGIN {{FS=\"\t\";OFS=\"\t\"}};{{print $2}}' | sort | uniq -c | awk 'BEGIN {{OFS=\"\t\"}};{{print $1,$2}}' > {}".format(tags_file,index_file),shell=True)
    return index_file

def reverseComplement(fq,dest="default"):
    if(dest=="default"):
        temp = '.'.join(fq.split('.')[0:-1])+"_rc.fastq"
        dest = temp.split('/')[-1]
    InFastq = SeqIO.parse(fq,"fastq")
    OutFastq = open(dest,"w",0)
    for record in InFastq:
        rc_record = record.reverse_complement(id=record.id+"_rc")
        SeqIO.write(rc_record,OutFastq,"fastq")
    OutFastq.close()
    return dest

def TrimAdapter(fq,adapter,trimfq="default",bar_len=12,tags_file="default",trim_err="default",start_trim=1,keepFailed=True):
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    if(trim_err=="default"):
        temp = '.'.join(fq.split('.')[0:-1]) + '.err.fastq'
        trim_err = temp.split('/')[-1]
    if(trimfq == "default"):
        temp = '.'.join(fq.split('.')[0:-1]) + '.trim.fastq'
        trimfq = temp.split('/')[-1]
    if(tags_file == "default"):
        temp = '.'.join(fq.split('.')[0:-1]) + '.tags.fastq'
        tags_file = temp.split('/')[-1]
    errOpen = open(trim_err,"w",0)
    tagsOpen = open(tags_file,"w",0)
    trimOpen = open(trimfq,"w",0)
    InFastq = SeqIO.parse(fq,"fastq")
    AdaptLen=len(adapter)
    TotalTrim=AdaptLen+bar_len+start_trim
    print("Adapter Length is {}".format(AdaptLen))
    for record in InFastq:
        pre_tag = SeqRecord(
                Seq(str(record.seq)[0:bar_len],"fastq"), \
                id=record.id)
        pre_tag.letter_annotations['phred_quality']=record.letter_annotations['phred_quality'][0:bar_len]
        '''
        if adapter not in pre_tag.seq:
            print("I'm sorry, but your adapter sequence is not in the tag. I will write this to an error fastq, which you are free to use or discard at your discretion")
            SeqIO.write(record,errOpen,"fastq")
            continue
        '''
        SeqIO.write(pre_tag,tagsOpen,"fastq")
        post_tag = SeqRecord(
                Seq(str(record.seq)[TotalTrim:],"fastq"), id=record.id,description=record.description+" ###" + str(record.seq[0:bar_len])) 
        post_tag.letter_annotations['phred_quality']=record.letter_annotations['phred_quality'][TotalTrim:]
        SeqIO.write(post_tag,trimOpen,"fastq")
    tagsOpen.close()
    trimOpen.close()
    errOpen.close()
    return(tags_file,trimfq)

def AdapterLoc(fq,adapter,bar_len=12,keepFailed=True):
    InFastq=SeqIO.parse(fq,"fastq")
    Tpref = '.'.join(fq.split('.')[0:-1])
    Prefix=Tpref.split('/')[-1]
    StdFilename, ElseFilename, ElseLoc = Prefix+'.{}.fastq'.format("adapter"+str(bar_len)), Prefix+ '.else.fastq', Prefix + '.else.supp'
    StdFastq=open(StdFilename,'w',0) #Adapter at expected Location
    ElseFastq=open(ElseFilename,'w',0) #Adapter sequence found elsewhere, even if it's simply part of the read itself
    ElseLocations=open(ElseLoc,'w',0)
    for read in InFastq:
        seq = str(read.seq)
        if(seq.find(adapter)==-1):
            if(keepFailed==False):
                raise NameError("Sanity check failure. Adapter Sequence Not Found. HTML 404")
            else:
                read.description+=" ###AdapterFail"
                SeqIO.write(read,StdFastq,"fastq")
        elif(seq[bar_len:bar_len+len(adapter)] == adapter):
            read.description+=" ###AdapterPass"
            SeqIO.write(read,StdFastq,"fastq") #Checks the expected adapter location. Avoiding find command, as sometimes the adapter sequence occurs before and at the expected location
        else:
            if(keepFailed==False):
                SeqIO.write(read,ElseFastq,"fastq")
                ElseLocations.write(repr(seq.find(adapter)) + "\t" + read.name + "\n")
            else:
                read.description=" ###AdapterFail"
                SeqIO.write(read,StdFastq,"fastq")
                ElseLocations.write(repr(seq.find(adapter)) + "\t" + read.name + "\n")
    StdFastq.close()
    ElseFastq.close()
    ElseLocations.close()
    return(StdFilename,ElseFilename)

def FastqRegex(fq,string,matchFile="default",missFile="default"):
    from subprocess import call
    if(matchFile=="default"):
        matchFile = ('.'.join(fq.split('.')[0:-1])+'.match.fastq').split('/')[-1]
    if(missFile=="default"):
        missFile = ('.'.join(fq.split('.')[0:-1])+'.miss.fastq').split('/')[-1]
    CommandStr = "cat {} | paste - - - - | grep '{}' | tr '\t' '\n' > {}".format(fq,string,matchFile)
    call(CommandStr,shell=True)
    CommandStr2 = "cat {} | paste - - - - | grep -v '{}' | tr '\t' '\n' > {}".format(fq,string,missFile)
    call(CommandStr2,shell=True)
    return(CommandStr,CommandStr2,matchFile,missFile)

def hamming(str1, str2):
    import operator
    from iterator import imap
    assert len(str1) == len(str2)
    #ne = str.__ne__  ## this is surprisingly slow
    #ne = operator.ne
    return sum(imap(operator.ne, str1, str2))

def sam_sort(insam,outsam):
    #Skip header and funnel all reads into a temp file
    import random
    output=open(outsam,'w+')
    tmpname='temp{}.txt'.format(random.randint(0,200))
    tmp=open(tmpname,'w',0)
    import subprocess
    command_str=str('grep -v "@SQ\|@PG\|VN:\|@HD" {}'.format(insam))
    print(command_str)
    subprocess.call(command_str,stdout=tmp,shell=True)
    tmp.close()
    #Save the header to the outsam
    command_header = 'grep "@SQ\|@PG\|@HD" {}'.format(insam)
    subprocess.call(command_header, stdout=output, shell=True)
    #sort the reads by query name
    tmp=open(tmpname,'r')
    command_str1=str('sort -k1,1 -t " " {}'.format(tmpname))
    print(command_str1)
    subprocess.call(command_str1,stdout=output,shell=True)
    output.close()
    tmp.close()
    subprocess.call('rm {}'.format(tmpname),shell=True)
    both_cmds=command_str+"\n"+command_str1
    return(both_cmds)

def Sam2Bam(insam,outbam):
    from subprocess import call
    output = open(outbam,'w',0)
    command_str='samtools view -Sbh {}'.format(insam,shell=True)
    print(command_str)
    call(command_str,stdout=output,shell=True)
    return(command_str,outbam)

def Bam2Sam(inbam,outsam):
    from subprocess import call
    output = open(outsam,'w',0)
    command_str='samtools view -h {}'.format(inbam)
    print(command_str)
    call(command_str,stdout=output,shell=True)
    return(command_str,outsam)

def fastx_trim(infq, outfq, n):
    import subprocess
    command_str = ['fastx_trimmer','-l',str(n),'-i',infq,'-o',outfq]
    print(command_str)
    subprocess.call(command_str)
    return(command_str)

def align_bwa(R1,R2,ref,opts,outsam):
    import subprocess
    opt_concat = ""
    if(opts== ""):
        opts='-t 4 -T 0 -v 1'
    output = open(outsam,'w',0)
    for i, opt_it in enumerate(opts.split()):
        opt_concat+=opt_it+" "
    command_str = 'bwa mem {} {} {} {}'.format(opt_concat,ref,R1,R2)
    #command_list = command_str.split(' ')
    print(command_str)
    subprocess.call(command_str, stdout=output,shell=True)
    output.close()
    return command_str;

def align_snap(R1,R2,ref,opts,outbam):
    import subprocess
    opt_concat = ""
    for i, opt_it in enumerate(opts.split()):
        opt_concat+=opt_it+" "
    command_str = "snap paired {} {} {} -o {} {}".format(ref,R1,R2,outbam,opt_concat)
    print(command_str)
    subprocess.call(command_str)
    return(command_str)

def align_bowtie2(R1,R2,ref,opts,outsam):
    import subprocess
    opt_concat=""
    output=open(outsam,'w',0) #Zero buffer size means it gets sent immediately straight to the file
    if(opts==""):
        opts='--threads 4 '
    if('--reorder' not in opts):
        opts+='--reorder '
    if('--mm' not in opts):
        opts+=' --mm '
    for i, opt_it in enumerate(opts.split()):
        opt_concat+=opt_it+" "
    command_str = 'bowtie2 {} --local --very-sensitive-local -x {} -1 {} -2 {}'.format(opt_concat,ref,R1,R2)
    print(command_str)
    #command_list=command_str.split(' ')
    subprocess.call(command_str, stdout=output,shell=True)
    output.close()
    return(command_str)

def fastq_sort(in_fastq,out_fastq):
    import subprocess
    outfile=open(out_fastq,'w')
    command_str='cat {} | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n"'.format(in_fastq)
    subprocess.call(command_str,stdout=outfile,shell=True)
    outfile.close()
    return(command_str)

if(__name__=="__main__"):
    main()
