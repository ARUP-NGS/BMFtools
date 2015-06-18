from __future__ import division
import numpy as np
from utilBMF.HTSUtils import pFastqProxy, pFastqFile, getBS
from itertools import groupby
from utilBMF.ErrorHandling import ThisIsMadness
from MawCluster.BCFastq import GetDescTagValue
import argparse

# create a generation from the pysam alignment object
# cimport pysam.calignmentfile
# ctypedef pysam.calignmentfile.AlignedSegment AlignedSegment_t
# inHandleCons = pysam.AlignmentFile(PathTowhatever)
# cons = inHandleCons.next()

parser=argparse.ArgumentParser(description="blah")
parser.add_argument("reads")
parser.add_argument("cons")
args=parser.parse_args()

consInFq = pFastqFile(args.cons)
readsInFq = pFastqFile(args.reads)

def createTransMatrix():
    transMatrix=[[0]*4,[0]*4]


nucIndex = {"A":0,"C":1,"T":3,"G":4,"N":5}
k=4
kmerDict={}
for bc4fq, fqRecGen, in groupby(readsInFq, key=getBS):
    consRead = consInFq.next()
    famCount = int(GetDescTagValue(consRead.comment,"FM"))
    if famCount <= 100:
        continue
    consBS = getBS(consRead)
    if(consBS != bc4fq):
        raise ThisIsMadness("Barcode mismatch, are both read files sorted?")
    pFqPrxList = list(fqRecGen)
    for rec in pFqPrxList:
        if(rec.sequence == consRead.sequence):
            continue
        for baseIndex in xrange(len(rec.sequence)-k):
            recKmer = rec.sequence[baseIndex:baseIndex+k]
            consKmer = consRead.sequence[baseIndex:baseIndex+k]
            recBase = rec.sequence[baseIndex+k]
            consBase = consRead.sequence[baseIndex+k]
            if(recKmer != consKmer):
                continue
            if recKmer not in kmerDict.keys():
                ()
                kmerDict[recKmer] = np.zeros((2,4,4),dtype=np.int)
            if recBase == consBase:
                kmerDict[recKmer][0][nucIndex[consBase]][nucIndex[recBase]]+=1
            else:
                kmerDict[recKmer][1][nucIndex[consBase]][nucIndex[recBase]]+=1

"""
            if recKmer not in kmerDict.keys():
                kmerDict[recKmer]={"A":[0,0],"C":[0,0],"T":[0,0],"G":[0,0],
                                   "N":[0,0],"total":[0,0]}
            if recBase == consBase:
                kmerDict[recKmer][recBase][0]+=1
                kmerDict[recKmer]["total"][0]+=1
            else:
                kmerDict[recKmer][recBase][1]+=1
                kmerDict[recKmer]["total"][1]+=1

nucList = ["total","A","T","C","G","N"]
out=open("Kmer_Bias.txt",'w')
for kmer in kmerDict:
    o = []
    o.append(kmer)
    #o.append(tr)
    #o.append(tw)
    for base in nucList:
        r = kmerDict[kmer][base][0]
        w = kmerDict[kmer][base][1]
        #o.append(r)
        #o.append(w)
        if w == 0:
            o.append(0)
        else:
            o.append(w/(r+w))
    o = [str(i) for i in o]
    o = "\t".join(o)+"\n"
    out.write(o)
out.close()
"""
