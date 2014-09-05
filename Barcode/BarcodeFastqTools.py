from Bio import SeqIO
import pysam
import logging

def BarcodeSort(inFastq,outFastq="default"):
    logging.info("Sorting fastq by barcode sequence.")
    from subprocess import call
    if(outFastq=="default"):
        outFastq = '.'.join(inFastq.split('.')[0:-1])+'.BS.fastq'
    BSstring = "cat {} | paste - - - - | grep 'HomingPass' | sed 's:###:###\t:g' |  awk 'BEGIN {{FS=OFS=\"\t\"}};{{print $3,$0}}' | sort -k1,1 | awk 'BEGIN {{OFS=FS=\"\t\"}};{{print $2,$3,$4,$5,$6,$7,$8,$9,$10}}' | sed 's:###\t:###:g' | sed 's:\t$::g' | sed 's:\t$::g' | tr '\t' '\n' > {}".format(inFastq,outFastq)
    call(BSstring,shell=True)
    logging.info("Command: {}".format(BSstring.replace("\t","\\t")))
    return outFastq

def compareFastqRecords(RecordList,stringency=0.9):
    from Bio.SeqRecord import SeqRecord
    seqs = [str(record.seq) for record in RecordList]
    max = 0
    Success = False
    for seq in seqs:
        #print("Seq: {}".format(str(seq)))
        numEq = sum(str(seq) == str(seqItem) for seqItem in seqs)
        if(numEq > max):
            max = numEq
            finalSeq = str(seq)
    frac = numEq * 1.0 / len(RecordList)
    #print("Fraction {}. Stringency: {}. Pass? {}.".format(frac,stringency,(frac>stringency)))
    if(frac > stringency):
        Success = True
    consolidatedRecord = SeqRecord(seq=finalSeq,id=RecordList[0].id,letter_annotations=RecordList[0].letter_annotations,name=RecordList[0].name,description=RecordList[0].description)
    return consolidatedRecord, Success

def HomingSeqLoc(fq,homing,bar_len=12):
    InFastq=SeqIO.parse(fq,"fastq")
    Tpref = '.'.join(fq.split('.')[0:-1])
    Prefix=Tpref.split('/')[-1]
    StdFilename, ElseFilename, ElseLoc = Prefix+'.{}.fastq'.format("homing"+str(bar_len)), Prefix+ '.else.fastq', Prefix + '.else.supp'
    StdFastq=open(StdFilename,'w',0) #Homing at expected Location
    ElseFastq=open(ElseFilename,'w',0) #Homing sequence found elsewhere, even if it's simply part of the read itself
    ElseLocations=open(ElseLoc,'w',0)
    for read in InFastq:
        seq = str(read.seq)
        if(seq.find(homing)==-1):
            read.description+=" ###HomingFail"
            SeqIO.write(read,StdFastq,"fastq")
        elif(seq[bar_len:bar_len+len(homing)] == homing):
            read.description+=" ###HomingPass"
            SeqIO.write(read,StdFastq,"fastq") #Checks the expected homing location. Avoiding find command, as sometimes the homing sequence occurs before and at the expected location
        else:
            read.description=" ###HomingFail"
            SeqIO.write(read,StdFastq,"fastq")
            ElseLocations.write(repr(seq.find(homing)) + "\t" + read.name + "\n")
    StdFastq.close()
    ElseFastq.close()
    ElseLocations.close()
    return StdFilename,ElseFilename

def fastq_sort(in_fastq,out_fastq):
    import subprocess
    outfile=open(out_fastq,'w')
    command_str='cat {} | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n"'.format(in_fastq)
    subprocess.call(command_str,stdout=outfile,shell=True)
    outfile.close()
    return(command_str)

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

def fastx_trim(infq, outfq, n):
    import subprocess
    command_str = ['fastx_trimmer','-l',str(n),'-i',infq,'-o',outfq]
    logging.info(command_str)
    subprocess.call(command_str)
    return(command_str)

def findProperPairs(infq1,infq2,index1="default",index2="default",outfq1="default",outfq2="default",outfqSingle="default"):
    logging.info("Now attempting to parse out proper pairs from reads which need to be aligned as single-end.")
    from Bio import SeqIO
    if(index1=="default"):
        raise ValueError("An index for the barcodes in reads1 set must be present.")
    if(index2=="default"):
        raise ValueError("An index for the barcodes in reads2 set must be present.")
    if(outfq1=="default"):
        outfq1 = '.'.join(infq1.split('.')[0:-1])+'.proper.fastq'
    if(outfq2=="default"):
        outfq2 = '.'.join(infq2.split('.')[0:-1])+'.proper.fastq'
    if(outfqSingle=="default"):
        outfqSingle = '.'.join(infq1.split('.')[0:-1])+'.solo.fastq'
    outfq1Handle = open(outfq1,'w')
    outfq2Handle = open(outfq2,'w')
    outfqSingleHandle = open(outfqSingle,'w')
    indexHandle1 = open(index1,"r")
    infq1Handle = open(infq1,"r")
    infq2Handle = open(infq2,"r")
    dictEntries = [line.split() for line in indexHandle1]
    BarDict1 = {}
    for entry in dictEntries:
        BarDict1[entry[1]]=entry[0]
    indexHandle2 = open(index2,"r")
    for read in infq2Handle:
        queryBC = read.description.split('###')[-2]
        try:
            temp = BarDict1[queryBC]
            SeqIO.write(read, outfq2Handle, "fastq")
        except KeyError:
            SeqIO.write(read,outfqSingleHandle,"fastq")
             
    dictEntries = [line.split() for line in indexHandle2]
    BarDict2 = {}
    for entry in dictEntries:
        BarDict2[entry[1]]=entry[0]
    for read in infq1Handle:
        queryBC = read.description.split('###')[-2]
        try:
            temp = BarDict2[queryBC]
            SeqIO.write(read, outfq1Handle, "fastq")
        except KeyError:
            SeqIO.write(read,outfqSingleHandle,"fastq")
    outfq1Handle.close()
    outfq2Handle.close()
    infq1Handle.close()
    infq2Handle.close()
    indexHandle1.close()
    indexHandle2.close()
    
    return outfq1,outfq1,outfqSingle

def GenerateSingleBarcodeIndex(tags_file,index_file="default"):
    from subprocess import call
    if(index_file=="default"):
        index_file = '.'.join(tags_file.split('.')[0:-1]) + ".barIdx"
    commandStr ="cat {} | sed 's:###: ###:g' | paste - - - - | grep -v \"HomingFail\" | awk 'BEGIN {{FS=\"\t\";OFS=\"\t\"}};{{print $2}}' | sort | uniq -c | awk 'BEGIN {{OFS=\"\t\"}};{{print $1,$2}}' > {}".format(tags_file,index_file)
    logging.info("CommandStr = {}".format(commandStr.replace("\t", "\\t")))
    call(commandStr,shell=True)
    return index_file


def GetFamilySizeSingle(trimfq,BarcodeIndex,outfq="default",singlefq="default"):
    logging.info("Getting family sizes for all of the grouped families.")
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
        try:
            famSize = BarDict[readTag]
        except KeyError:
            famSize= 0
        newRead.description = read.description+" ###"+str(famSize)
        #print("famSize = _{}_".format(str(famSize)))
        #print("The value of this comparison to 1 is {}".format(str(famSize=="1")))
        if(str(famSize) == "0"):
            continue
        if(str(famSize) == "1"):
            SeqIO.write(newRead,singlefqBuffer,"fastq")
            SeqIO.write(newRead,outfqBuffers,"fastq")
        else:
            ReadsWithFamilies+=1
            SeqIO.write(newRead,outfqBuffers,"fastq")
    return outfq,TotalReads, ReadsWithFamilies

def mergeSequencesFastq(fq1,fq2,output="default"):
    if(output=="default"):
        output = fq1.split('.')[0] + '.merged.fastq'
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    reads1=SeqIO.parse(fq1,"fastq")
    reads2=SeqIO.parse(fq2,"fastq")
    outFastq = open(output,"w",0)
    while True:
        try:
            read1 = reads1.next()
            read2 = reads2.next()
            outread = read1
            outread=SeqRecord(
                Seq(str(read1.seq)+str(read2.seq),"fastq"), id=read1.id,description=read1.description) 
            outread.letter_annotations['phred_quality']=read1.letter_annotations['phred_quality'] + read2.letter_annotations['phred_quality'] 
            SeqIO.write(outread, outFastq, "fastq")
        except StopIteration:
            break
    reads1.close()
    reads2.close()
    outFastq.close()
    return output

def pairedFastqConsolidate(fq1,fq2,outFqPair1="default",outFqPair2="default",stringency=0.9):
    logging.info("Now consolidating paired-end reads.")
    if(outFqPair1=="default"):
        outFqPair1= '.'.join(fq1.split('.')[0:-1]) + 'cons.fastq'
    if(outFqPair2=="default"):
        outFqPair2= '.'.join(fq2.split('.')[0:-1]) + 'cons.fastq'
    inFq1 = SeqIO.parse(fq1,'fastq')
    inFq2 = SeqIO.parse(fq2,'fastq')
    outputHandle1 = open(outFqPair1,'w')
    outputHandle2 = open(outFqPair2,'w')
    workingBarcode = ""
    workingSet1 = []
    workingSet2 = []
    for fqRec in inFq1:
        fqRec2 = inFq2.next()
        barcode4fq1 = fqRec.description.split("###")[-2].strip()
        #print("Working barcode: {}. Current barcode: {}.".format(workingBarcode,barcodeRecord))
        #print("name of read with this barcode: {}".format(record.qname))
        #print("Working set: {}".format(workingSet1))
        if("AAAAAAAAAA" in barcode4fq1 or "TTTTTTTTTT" in barcode4fq1 or "CCCCCCCCCC" in barcode4fq1 or "GGGGGGGGGG" in barcode4fq1):
            continue
        if(int(fqRec.description.split('###')[-1].strip()) < 2):
            continue
        if(workingBarcode == ""):
            workingBarcode = barcode4fq1
            workingSet1 = []
            workingSet1.append(fqRec)
            workingSet2 = []
            workingSet2.append(fqRec2)
            continue
        elif(workingBarcode == barcode4fq1):
            workingSet1.append(fqRec)
            workingSet2.append(fqRec2)
            continue
        elif(workingBarcode != barcode4fq1):
            mergedRecord1, success1 = compareFastqRecords(workingSet1,stringency=float(stringency))
            mergedRecord2, success2 = compareFastqRecords(workingSet2,stringency=float(stringency))
            if(success1 == True):
                SeqIO.write(mergedRecord1, outputHandle1, "fastq")
            if(success2 == True):
                SeqIO.write(mergedRecord2, outputHandle2, "fastq")
            workingSet1 = []
            workingSet1.append(fqRec)
            workingSet2 = []
            workingSet2.append(fqRec2)
            workingBarcode = barcode4fq1
            continue
        else:
            raise RuntimeError("No idea what's going on. This code should be unreachable")
    inFq1.close()
    inFq2.close()
    outputHandle1.close()
    outputHandle2.close()
    return outFqPair1,outFqPair2

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

def singleFastqConsolidate(fq,outFq="default",stringency=0.9):
    if(outFq=="default"):
        outFq= '.'.join(fq.split('.')[0:-1]) + 'cons.fastq'
    inFq = SeqIO.parse(fq,'fastq')
    outputHandle = open(outFq,'w')
    workingBarcode = ""
    workingSet = []
    for fqRec in inFq:
        barcode4fq = fqRec.description.split("###")[-2].strip()
        #print("Working barcode: {}. Current barcode: {}.".format(workingBarcode,barcodeRecord))
        #print("name of read with this barcode: {}".format(record.qname))
        #print("Working set: {}".format(workingSet))
        if("AAAAAAAAAA" in barcode4fq or "TTTTTTTTTT" in barcode4fq or "CCCCCCCCCC" in barcode4fq or "GGGGGGGGGG" in barcode4fq):
            continue
        if(int(fqRec.description.split('###')[-1].strip()) < 2):
            continue
        if(workingBarcode == ""):
            workingBarcode = barcode4fq
            workingSet = []
            workingSet.append(fqRec)
            continue
        elif(workingBarcode == barcode4fq):
            workingSet.append(fqRec)
            continue
        elif(workingBarcode != barcode4fq):
            mergedRecord, success = compareFastqRecords(workingSet,stringency=float(stringency))
            if(success == True):
                SeqIO.write(mergedRecord, outputHandle, "fastq")
            workingSet = []
            workingSet.append(fqRec)
            workingBarcode = barcode4fq
            continue
        else:
            raise RuntimeError("No idea what's going on. This code should be unreachable")
    inFq.close()
    outputHandle.close()
    return outFq

def TrimHoming(fq,homing,trimfq="default",bar_len=12,tags_file="default",trim_err="default",start_trim=1):
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
    tagsOpen = open(tags_file,"w",0)
    trimOpen = open(trimfq,"w",0)
    InFastq = SeqIO.parse(fq,"fastq")
    HomingLen=len(homing)
    TotalTrim=HomingLen+bar_len+start_trim
    logging.info("Homing Length is {}".format(HomingLen))
    for record in InFastq:
        pre_tag = SeqRecord(
                Seq(str(record.seq)[0:bar_len],"fastq"), \
                id=record.id,description=record.description)
        pre_tag.letter_annotations['phred_quality']=record.letter_annotations['phred_quality'][0:bar_len]
        '''
        if homing not in pre_tag.seq:
            logging.info("I'm sorry, but your homing sequence is not in the tag. I will write this to an error fastq, which you are free to use or discard at your discretion")
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
    return(tags_file,trimfq)