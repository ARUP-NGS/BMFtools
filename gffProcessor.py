#!/mounts/anaconda/bin/python
import argparse
def main():

    parser=argparse.ArgumentParser()
    parser.add_argument('-i','--gff', help="Provide your gff file", metavar = ('Input GFF'),required=True)
    parser.add_argument('-o','--outfile', help="Set output file for processed GFF. Default = stdout", nargs = 1, metavar = ('Output GFF'))
    args=parser.parse_args()
    #GFF,p1,p2=FetchGFF(args.gff)
    GFF=FetchGFF(args.gff)
    print("GFF loaded. Now separating")
    p1,p2=[],[]
    for rec in GFF:
        p1.append(GFF[:8])
        p2.append(GFF[8:])
    del GFF
    print("GFF separated, first object cleared. Now modifying")
    outfile=args.outfile
    NewP1=ExpandGffP1(p1,2000)
    del p1
    print("Now writing to file {}".format(outfile))
    catGff(NewP1,p2,outfile)
    return


def FetchGFF(gff): 
    rawContent = [line.strip() for line in open(gff)]
    tsv,ssv,content=[],[],[] #Tab-separated,semicolon-separated,full-record
    for count, rec in enumerate(rawContent):
        if(rec[0]=="#"):
            print("Line # {} begins with a comment. Skipping line.".format(count))
            continue
        else:
            recSplit=rec.split("\t")
            temp_ssv=recSplit[-1].split(";") #Snag the last entry for the semicolon-separated data
            temp_tsv=recSplit[0:-1]
            tsv.append(temp_tsv)
            ssv.append(temp_ssv)
            for entry in temp_ssv:
                temp_tsv.append(entry)
            content.append(temp_tsv) #TODO: clean this up. Apparently, tsv and content are redundant at the end of execution
    '''
    p1=[]
    p2=[]
    for rec in content:
        p1.append(rec[:8])
        p2.append(rec[8:])
    return content,p1,p2
    '''
    return content #TODO: try expanding content to splitting it up.
#class GffRecord:

#Prints out a proper GFF from the representation of a GFF record in this program
def catGff(p1,p2,outfile):
    handle=open(outfile,"w",0)
    for rec1,rec2 in zip(p1,p2):
        LineToWrite = "\t".join(rec1) + "\t" + ";".join(rec2) + "\n"
        handle.write(LineToWrite)
    return

def ExpandGffP1(p1,len):
    NewP1 = []
    for rec in p1:
        gff_hold=rec
        print("gff_hold 4 is {}".format(gff_hold[4]))
        print("rec is {}".format(rec))
        gff_hold[4] = str(int(gff_hold[4])-len)
        gff_hold[5] = str(int(gff_hold[5])-len)
        NewP1.append(gff_hold)
    return NewP1

if(__name__=="__main__"):
    main()

