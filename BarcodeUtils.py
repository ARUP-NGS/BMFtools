#!/mounts/anaconda/bin/python

from Bio import SeqIO
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--fq', help="Provide your fastq file", nargs = 1, metavar = ('reads'),required=True)
    parser.add_argument('-b','--adapter', help="Adapter for samples. If not set, defaults to CAGT.",metavar=('Adapter'),default="CAGT")
    parser.add_argument('-l','--bar_len', help="Length of anticipated barcode. Defaults to 12.", default=12, type=int)
    parser.add_argument('-p','--pair', help="Unpaired (0), Paired End 1 (1), Paired End 2 (2). Defaults to 0.", default=0, type=int)
    args=parser.parse_args()
    adapter=args.adapter
    fq = args.fq[0]
    Regex1,Regex2,Hits,Misses=FastqRegex(fq,adapter)
    print("Regex operation complete. Now locating adapter sequence.")
    StdFilenames,ElseFilenames=AdapterLoc(Hits,adapter=adapter)
    print("Adapter sequences located. Hits and misses (correct location vs. incorrect location or not found) parsed out.")
    print("Now removing the adapter and the barcode.")
    tags, trimfq = TrimAdapter(StdFilenames,adapter)
    BarcodeIndex = GenerateBarcodeIndex(tags)
    if(args.pair == 2):
        revCompTrim = reverseComplement(trimfq)
        print("revCompTrim location: " + revCompTrim)

        #Now, "tags" is the file location for the tags. I will use it to make a table of the counts for each family.
    return

def GetFamilySize(trimfq,BarcodeIndex):


def GenerateBarcodeIndex(tags_file,index_file="default"):
    from subprocess import call
    if(index_file=="default"):
        index_file = '.'.join(tags_file.split('.')[0:-1]) + ".barIdx"
    call("cat {} | paste - - - - | awk 'BEGIN {FS=\"\t\";OFS=\"\t\"};{print $2}' | sort | uniq -c | awk 'BEGIN {OFS=\"\t\"};{print $1,$2}' > {}".format(tags_file,index_file),shell=True)
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

def TrimAdapter(fq,adapter,trimfq="default",bar_len=12,tags_file="default",trim_err="default",start_trim=1):
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
        if adapter not in pre_tag.seq:
            print("I'm sorry, but your adapter sequence is not in the tag. I will write this to an error fastq, which you are free to use or discard at your discretion")
            SeqIO.write(record,errOpen,"fastq")
            continue
        SeqIO.write(pre_tag,tagsOpen,"fastq")
        post_tag = SeqRecord(
                Seq(str(record.seq)[TotalTrim:],"fastq"), id=record.id + "\t###" + str(record.seq)[0:bar_len],description="") 
        post_tag.letter_annotations['phred_quality']=record.letter_annotations['phred_quality'][TotalTrim:]
        SeqIO.write(post_tag,trimOpen,"fastq")
    tagsOpen.close()
    trimOpen.close()
    errOpen.close()
    return(tags_file,trimfq)

def AdapterLoc(fq,adapter,bar_len=12):
    InFastq=SeqIO.parse(fq,"fastq")
    Tpref = '.'.join(fq.split('.')[0:-1])
    Prefix=Tpref.split('/')[-1]
    StdFilename, ElseFilename, ElseLoc = Prefix+'.{}.fastq'.format(bar_len), Prefix+ '.else.fastq', Prefix + '.else.supp'
    StdFastq=open(StdFilename,'w',0) #Adapter at expected Location
    ElseFastq=open(ElseFilename,'w',0) #Adapter sequence found elsewhere, even if it's simply part of the read itself
    ElseLocations=open(ElseLoc,'w',0)
    for read in InFastq:
        seq = str(read.seq)
        if(seq.find(adapter)==-1):
            raise NameError("Sanity check failure. Adapter Sequence Not Found. HTML 404")
        elif(seq.find(adapter)==bar_len):
            SeqIO.write(read,StdFastq,"fastq")
        elif(seq[bar_len:bar_len+len(adapter)] == adapter):
            SeqIO.write(read,StdFastq,"fastq") #Sometimes there is a match for the adapter before the adapter occurs. Check the reads which don't find the adapter at the expected location directly
        else:
            SeqIO.write(read,ElseFastq,"fastq")
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

'''
def SumStrList(list):
    a = ""
    for i in list:
        a+=i
    return a
'''

if(__name__=="__main__"):
    main()
