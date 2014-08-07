#!/mounts/anaconda/bin/python
import argparse
#This program is motivated by the problem of deducing the effects of PCR duplication
#and its error introduction on 


def main():

    parser=argparse.ArgumentParser()
    parser.add_argument("-b","--batched",help="Calculate frequencies for a range of numbers of replicates in a duplicate family.", default=False)
    parser.add_argument("--lb",help="Lower bound for the range of family sizes. Only used for batch. Defaults to 5", default=5)
    parser.add_argument("--ub",help="Upper bound for the range of family sizes. Only used for batch. Defaults to 100", default=100)
    parser.add_argument("--ss",help="Step size for the range of family sizes.", default="5")
    parser.add_argument("-i","--iterations",help="Number of Monte Carlo iterations. Used for both batch and batch==False",default=100000,type=int)
    parser.add_argument("-n","--size-duplicate-family",help="Number of duplications in the duplicate family, one less than the number in the family. Only used for batch==False",default=10, type=int)
    parser.add_argument("-o","--outfile",help="File to write output.",default="")
    args=parser.parse_args()
    size=args.size_duplicate_family
    outfile=args.outfile

    if(args.batch==False):
        finalSumHist,outfile = filesToUse(size_duplicate_family=size,iterations=args.iterations,outfile=outfile)
    if(args.batch==True):
        values = 
    

    return

def filesToUse(size_duplicate_family=10,iterations=100000,random=False,outfile="defaultOutfile.omgz"):
    '''
    parser=argparse.ArgumentParser()
    parser.add_argument("-n","--size-duplicate-family",help="Number of duplications in the duplicate family, one less than the number in the family",default=10, type=int)
    parser.add_argument("-i","--iterations",help="Number of Monte Carlo iterations.",default=100000,type=int)
    parser.add_argument("-r","--random",help="Randomly select the number of reads for a duplicate family",default=False,type=bool)
    parser.add_argument("-o","--outfile",help="File to write output.",default="")
    args=parser.parse_args()
    size=int(args.size_duplicate_family)
    outfile=args.outfile
    '''
    size=size_duplicate_family
    if(outfile == ""):
        outfile = "PCR_Duplication_Simulation_{}_InFamily_{}_Iterations.txt".format(size,iterations)
    countMC=0;
    sumHist=[0] * size
    print("Beginning loop")
    while(countMC < int(iterations)):
        iterationArray = dupDistribution(size)
        for count,entry in enumerate(iterationArray):
            sumHist[count]+=entry
        countMC+=1
    outFile=open(outfile,"w",0)
    outFile.write("@Description: PCR Duplication Simulation Results.\n")
    outFile.write("@Header-Copied_#_times\tNumberOfReadsMatchingSaidDescription\n")
    for count,dup in enumerate(sumHist):
        outFile.write("{}\t{}\n".format(count,dup))
    outFile.write(repr(sumHist))
    return sumHist, outFile


def dupDistribution(size):
    import random
    random.seed()
    dupHist=[0.]*size
    #dupHist[0],dupHist[1]=1,1
    dupHist[0]=1
    while(sum(dupHist)<size+1):
        cdf=[]
        for count, i in enumerate(dupHist):
            cdf.append(sum(dupHist[0:count]))
        sum_cdf=cdf[-1]
        bins=[(dupHist[x])/sum_cdf for x in range(size)]
        dropTheNeedle=random.random()
        for count, bin in enumerate(bins):
            if(dropTheNeedle < bin):
                dupHist[count+1]+=1
                break
    return dupHist  
            
    print("Completing distribution of duplication.")
    print(repr(dupHist))
    return dupHist 

if(__name__=="__main__"):
    main()

