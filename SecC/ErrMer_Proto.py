from __future__ import division
import numpy as np
from utilBMF.HTSUtils import pFastqProxy, pFastqFile, getBS
from itertools import groupby
from utilBMF.ErrorHandling import ThisIsMadness
from MawCluster.BCFastq import GetDescTagValue
import argparse
import operator
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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


nucIndex = {"A":0,"C":1,"T":2,"G":3}
nucList = ["A","C","T","G","N"]
k=4
kmerDict={}
kmerTotals={}
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
            if recBase == "N":
                continue
            if recKmer not in kmerDict.keys():
                kmerDict[recKmer] = np.zeros((4,4),dtype=np.float)
                kmerTotals[recKmer] = [0,0]
            kmerDict[recKmer][nucIndex[recBase]][nucIndex[consBase]]+=1.0
            kmerTotals[recKmer][1]+=1
            if recBase != consBase:
                kmerTotals[recKmer][0]+=1

wrongKmers = {kmer:kmerTotals[kmer][0]/kmerTotals[kmer][1] for kmer
              in kmerTotals}
sortedKmers = sorted(wrongKmers.items(), key=operator.itemgetter(1),
                     reverse=True)
heatedKmers = {kmer:kmerDict[kmer]/kmerDict[kmer].sum(axis=0)[:,None] for
               kmer in kmerDict}

for kmer, wrong in sortedKmers[:10]:
    fig, ax = plt.subplots()
    plt.suptitle(kmer,fontsize=24)
    data = heatedKmers[kmer]
    np.fill_diagonal(data,0)
    heatmap = ax.pcolor(data, cmap=plt.cm.Blues, edgecolors='k')
    ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
    ax.set_yticks(np.arange(data.shape[1])+0.5, minor=False)
    ax.set_xticklabels(nucList, minor=False)
    ax.set_yticklabels(nucList, minor=False)
    for y in range(data.shape[0]):
        for x in range(data.shape[1]):
            plt.text(x+0.5, y+0.5, str(heatedKmers[kmer][y, x])[:5],
                     horizontalalignment='center',
                     verticalalignment='center',
                     )
    plt.show()
