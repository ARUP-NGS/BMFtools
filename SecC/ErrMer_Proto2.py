from __future__ import division
import numpy as np
from utilBMF.HTSUtils import pFastqProxy, pFastqFile, getBS, RevCmp
from itertools import groupby
from utilBMF.ErrorHandling import ThisIsMadness
from MawCluster.BCFastq import GetDescTagValue
import argparse
import operator
import pysam
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
mcoptBS = operator.methodcaller("opt","BS")
# cimport pysam.calignmentfile
# ctypedef pysam.calignmentfile.AlignedSegment AlignedSegment_t
nucList = ["A", "T", "C", "G", "N"]


def getArgs():
    parser = argparse.ArgumentParser(description="blah")
    parser.add_argument("mode", help="current options: heatmap, cycleError")
    parser.add_argument("reads1")
    parser.add_argument("reads2")
    parser.add_argument("bam", help="aligned bam if calling heater, calMD bam"
                                     "if calling cycleError")
    parser.add_argument("kmerLen")
    return parser.parse_args()


def HeatMaps(heatedKmers, topKmers):
    for kmer, wrong in topKmers:
        fig, ax = plt.subplots()
        plt.suptitle(kmer+"  "+str(wrong), fontsize=24)
        data = heatedKmers[kmer]
        # np.fill_diagonal(data,0)
        heatmap = ax.pcolor(data, cmap=plt.cm.Blues, edgecolors='k')
        ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
        ax.set_yticks(np.arange(data.shape[1])+0.5, minor=False)
        ax.set_xticklabels(nucList, minor=False)
        ax.set_yticklabels(nucList, minor=False)
        plt.xlabel("Consensus Base")
        plt.ylabel("Record Base")
        for y in range(data.shape[0]):
            for x in range(data.shape[1]):
                plt.text(x+0.5, y+0.5, str(heatedKmers[kmer][y, x])[:5],
                         horizontalalignment='center',
                         verticalalignment='center',)
        plt.savefig(kmer)


def recsConsCompare(recGen, consRead, kmerDict, kmerTotals, kmerQuals, k):
    """Compares read and consensus to build a dictionary of context specific
        errors.    SOON TO BE DEPRICATED"""
    nucIndex = {"A": 0, "C": 1, "T": 2, "G": 3}
    for rec in recGen:
        recSeq = rec.sequence
        recQual = rec.getQualArray()
        if consRead.is_reverse:
            conSeq = RevCmp(consRead.query_sequence)
            conSeq = consRead.query_qualities[::-1]
        else:
            conSeq = consRead.query_sequence
            conQual = consRead.query_qualities
        if(recSeq == conSeq):
            continue
        for baseIndex in xrange(len(rec.sequence)-k):
            recKmer = recSeq[baseIndex:baseIndex+k]
            consKmer = conSeq[baseIndex:baseIndex+k]
            if(recKmer != consKmer):
                continue
            recBase = recSeq[baseIndex+k]
            consBase = conSeq[baseIndex+k]
            qualScore = conQual
            if recBase == "N":
                continue
            try:
                kmerDict[recKmer]
                pass
            except KeyError:
                kmerDict[recKmer] = np.zeros((4, 4), dtype=np.float)
                kmerTotals[recKmer] = [0, 0]
                kmerQuals[recKmer] = np.zeros((4, 4), dtype=np.float)
            kmerDict[recKmer][nucIndex[consBase]][nucIndex[recBase]] += 1.0
            # kmerQuals[recKmer][nucIndex[consBase]][nucIndex[recBase]][0] += conQual
            # kmerQuals[recKmer][nucIndex[consBase]][nucIndex[recBase]][1] += 1.0
            kmerTotals[recKmer][1] += 1.0
            if recBase != consBase:
                kmerTotals[recKmer][0] += 1.0

def errorFinder(recGen, mdRead, readError, readObs):
    print "ding!"
    for rec in recGen:
        recSeq = rec.sequence
        if mdRead.is_reverse:
            mdSeq = RevCmp(mdRead.query_sequence)
        else:
            mdSeq = mdRead.query_sequence
        for baseIndex in xrange(len(recSeq)):
            print mdSeq[baseIndex]
            print recSeq[baseIndex]

def cycleError(args):
    print "ding!"
    readsInFq1 = pFastqFile(args.reads1)
    readsInFq2 = pFastqFile(args.reads2)
    MDbamGen = groupby(pysam.AlignmentFile(args.bam), key=mcoptBS)
    bgn = MDbamGen.next
    r2GenGen = groupby(readsInFq2, key=getBS)
    r2ggn = r2GenGen.next
    read1error = []
    read2error = []
    read1obs = []
    read2obs = []
    for bc4fq1, fqRecGen1 in groupby(readsInFq1, key=getBS):
        MDbs, mdReads = bgn()
        bc4fq2, fqRecGen2 = r2ggn()
        mdRead1 = None
        mdRead2 = None
        if(bc4fq1 != bc4fq2):
            raise ThisIsMadness("Fastq barcodes don't match")
        if(MDbs != bc4fq1):
            raise ThisIsMadness("Bam Barcode mismatch")
        for cRead in mdReads:
            if cRead.is_secondary:
                continue
            if cRead.is_supplementary:
                continue
            if cRead.is_read1:
                mdRead1 = cRead
            if cRead.is_read2:
                mdRead2 = cRead
        if mdRead1 is None or mdRead2 is None:
            print "Did not find both reads in bam"
            continue
        if mdRead1.opt("FM") <= 100:
            continue
        if mdRead1.mapping_quality <= 20 or mdRead2.mapping_quality <= 20:
            continue
        errorFinder(fqRecGen1, mdRead1, read1error, read1obs)
        errorFinder(fqRecGen2, mdRead2, read2error, read2obs)
    read1error = np.array(read1error)
    read1obs = np.array(read1obs)
    read2error = np.array(read2error)
    read2obs = np.array(read2obs)
    read1PropError = read1error/read1obs
    read2PropError = read2error/read2obs

def heater(args):
    readsInFq1 = pFastqFile(args.reads1)
    readsInFq2 = pFastqFile(args.reads2)
    bamGen = groupby(pysam.AlignmentFile(args.bam), key=getBS)
    k = int(args.kmerLen)
    bgn = bamGen.next
    r2GenGen = groupby(readsInFq2, key=getBS)
    r2ggn = r2GenGen.next
    kmerDict = {}
    kmerTotals = {}
    kmerQuals = {}
    for bc4fq1, fqRecGen1, in groupby(readsInFq1, key=getBS):
        consBS, consReads = bgn()
        bc4fq2, fqRecGen2 = r2ggn()
        consRead1 = None
        consRead2 = None
        if(bc4fq1 != bc4fq2):
            raise ThisIsMadness("Fastq barcodes don't match")
        if(consBS != bc4fq1):
            raise ThisIsMadness("Bam Barcode mismatch")
        for cRead in consReads:
            if cRead.is_secondary:
                continue
            if cRead.is_supplementary:
                continue
            if cRead.is_read1:
                consRead1 = cRead
            if cRead.is_read2:
                consRead2 = cRead
        if consRead1 is None or consRead2 is None:
            print "Did not find both reads in bam"
        if consRead1.opt("FM") <= 100:
            continue
        recsConsCompare(fqRecGen1, consRead1, kmerDict, kmerTotals, kmerQuals,
                        k)
        recsConsCompare(fqRecGen2, consRead2, kmerDict, kmerTotals, kmerQuals,
                        k)
    wrongKmers = {kmer: kmerTotals[kmer][0]/kmerTotals[kmer][1] for kmer
                  in kmerTotals}
    sortedKmers = sorted(wrongKmers.items(), key=operator.itemgetter(1),
                         reverse=True)
    heatedKmers = {kmer: kmerDict[kmer]/kmerDict[kmer].sum(axis=0) for kmer in
                   kmerDict}
    HeatMaps(heatedKmers, sortedKmers[:10])


def main():
    print "ding!"
    args = getArgs()
    if args.mode not in ["heatmap", "cycleError"]:
        raise ThisIsMadness("invalid mode")
    if args.mode == "heater":
        heater(args)
    if args.mode == "cycleError":
        cycleError(args)


if __name__ == "__main__":
    main()
