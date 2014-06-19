#!/mounts/anaconda/bin/python

from Bio import SeqIO
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--fq', help="Provide your fastq file", nargs = 1, metavar = ('reads'),required=True)
    parser.add_argument('-b','--adapter', help="Adapter for samples. If not set, defaults to CAGT.",metavar=('Adapter'),default="CAGT")
    parser.add_argument('-l','--bar_len', help="Length of anticipated barcode.", default=12)
    args=parser.parse_args()
    adapter=args.adapter
    fq = args.fq[0]
    Regex1,Regex2,Hits,Misses=FastqRegex(fq,adapter)
    print("Regex operation complete. Now locating adapter sequence.")
    StdFilenames,ElseFilenames=AdapterLoc(Hits,adapter=adapter)
    print("Adapter sequences located. Hits and misses (correct location vs. incorrect location or not found) parsed out.")
    print("Now removing the adapter and the barcode.")
    tags, trimfq = TrimAdapter(Hits,adapter)
    return

def TrimAdapter(fq,adapter,trimfq="default",bar_len=12,tags_file="default"):
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    if(trimfq == "default"):
        trimfq = '.'.join(fq.split('.')[0:-1]) + '.trim.fastq'
    if(tags_file == "default"):
        tags_file = '.'.join(fq.split('.')[0:-1]) + '.tags.fastq'
    tagsOpen = open(tags_file,"w",0)
    trimOpen = open(trimfq,"w",0)
    InFastq = SeqIO.parse(fq,"fastq")
    AdaptLen=len(adapter)
    TotalTrim=AdaptLen+bar_len
    print("Adapter Length is {}".format(AdaptLen))
    for record in InFastq:
        pre_tag = SeqRecord(
                Seq(str(record.seq)[0:TotalTrim],"fastq"), \
                id=record.id, name=record.name)
        pre_tag.letter_annotations['phred_quality']=record.letter_annotations['phred_quality'][0:TotalTrim]
        SeqIO.write(pre_tag,tagsOpen,"fastq")
        post_tag = SeqRecord(
                Seq(str(record.seq)[TotalTrim:],"fastq"), \
                id=record.id, name=record.name)
        post_tag.letter_annotations['phred_quality']=record.letter_annotations['phred_quality'][TotalTrim:]
        SeqIO.write(post_tag,trimOpen,"fastq")
    tagsOpen.close()
    trimOpen.close()
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
        if(seq.find(adapter)==bar_len):
            SeqIO.write(read,StdFastq,"fastq")
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
        matchFile = '.'.join(fq.split('.')[0:-1])+'.match.fastq'
    if(missFile=="default"):
        missFile = '.'.join(fq.split('.')[0:-1])+'.miss.fastq'
    CommandStr = "cat {} | paste - - - - | grep '{}' | tr '\t' '\n' > {}".format(fq,string,matchFile)
    call(CommandStr,shell=True)
    CommandStr2 = "cat {} | paste - - - - | grep -v '{}' | tr '\t' '\n' > {}".format(fq,string,missFile)
    call(CommandStr2,shell=True)
    return(CommandStr,CommandStr2,matchFile,missFile)

def SumStrList(list):
    a = ""
    for i in list:
        a+=i
    return a

if(__name__=="__main__"):
    main()
