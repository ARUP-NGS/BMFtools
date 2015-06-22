from __future__ import division
import numpy as np
from utilBMF.HTSUtils import pFastqProxy, pFastqFile, getBS
from utilBMF.HTSUtils cimport nucList
from itertools import groupby
from utilBMF.ErrorHandling import ThisIsMadness
from MawCluster.BCFastq import GetDescTagValue
import argparse
import operator
import pysam
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
cimport pysam.calignmentfile
ctypedef pysam.calignmentfile.AlignedSegment AlignedSegment_t
mcoptBS = operator.methodcaller("opt","BS")

def getArgs():
    parser=argparse.ArgumentParser(description="blah")
    parser.add_argument("reads")
    parser.add_argument("consBam")
    parser.add_argument("kmer")
    return parser.parse_args()

def HeatMaps(heatedKmers,sortedKmers):
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
                         verticalalignment='center',)
        plt.show()

def ParseMDZ(mdzStr,start):
    """
    Returns the base position of mismatches between the current read and
    the reference, given the MD Z string and the start position for the read.
    """
    lastval = None
    oper = []
    bases = ["A","T","C","G","N"]
    counts = [str(i) for i in range(0,10)]+["~"]
    for item in mdzStr:
        curval = item
        if item in bases:
            itemType = "base"
        if item in counts:
            itemType = "count"
        if curval == lastval:
            raise NotImplementedError("Oh no, I haven't finished this...")


def main():
    args = getArgs()
    consBam = pysam.AlignmentFile(args.consBam)
    readsInFq = pFastqFile(args.reads)
    bamGen = groupby(consBam,key=mcoptBS)
    nucIndex = {"A":0,"C":1,"T":2,"G":3}
    nucList = ["A","C","T","G","N"]
    k=4
    kmerDict={}
    kmerTotals={}
    for bc4fq, fqRecGen, in groupby(readsInFq, key=getBS):
        consRead = consBam.next()
        consBS = bc4fq
        famCount = consRead.opt("FM")
        if famCount <= 100:
            continue
        if(consBS != bc4fq):
            raise ThisIsMadness("Barcode mismatch, are both read files"
                                "sorted?")
        pFqPrxList = list(fqRecGen)
        for rec in pFqPrxList:
            if(rec.sequence == consRead.sequence):
                continue
            for baseIndex in xrange(len(rec.sequence)-k):
                recKmer = rec.sequence[baseIndex:baseIndex+k]
                consKmer = consRead.query_sequence[baseIndex:baseIndex+k]
                recBase = rec.sequence[baseIndex+k]
                consBase = consRead.query_sequence[baseIndex+k]
                if(recKmer != consKmer):
                    continue
                if recBase == "N":
                    continue
                '''
                Replaced an "x in y" call with a dict try/except.
                if recKmer not in kmerDict.keys():
                '''
                try:
                    kmerDict[recKmer]
                    pass
                except KeyError:
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
    HeatMaps()


if __name__ == "__main__":
    main()
