#/mounts/anaconda/bin/python

import argparse

def main():

    parser=argparse.ArgumentParser()
    parser.add_argument("--decrement-only","-d",help="Boolean value. If you have already processed your maf files with a bash script and simply need to change from 1-based to 0-based coordinates, select this.",type=bool,default=True)
    parser.add_argument("--input","-i",help="Input file, whether MAF or processed to a one-off bedfile.")
    parser.add_argument("--output","-o",help="Output file. Defaults to a variation on the input name.",default="default")
    args=parser.parse_args()
    infile = args.input
    outfile = args.output
    if(outfile=="default"):
        outfile = '.'.join(infile.split('.')[0:-1])+".corr.bed"
    inhandle = open(infile,"r")
    outhandle = open(outfile,"w",0)
    for line in inhandle:
        bedrow = line.split('\t')
        try:
            bedrow[1] = str(int(bedrow[1])-1)
            bedrow[2] = str(int(bedrow[2])-1)
            bedwrite = '\t'.join(bedrow)
            outhandle.write(bedwrite)
        except ValueError:
            continue
    inhandle.close()
    outhandle.close()
    return    

if(__name__=="__main__"):
    main()
