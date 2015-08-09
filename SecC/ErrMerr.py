from __future__ import division
import numpy as np
from utilBMF.HTSUtils import pFastqProxy, pFastqFile, getBS, RevCmp
from itertools import groupby
from utilBMF.ErrorHandling import ThisIsMadness
from MawCluster.BCFastq import GetDescTagValue
import argparse
import operator
import pysam
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
mcoptBS = operator.methodcaller("opt","BS")
# cimport pysam.calignmentfile
# ctypedef pysam.calignedsegment.AlignedSegment AlignedSegment_t
nucList = ["A", "T", "C", "G", "N"]


def getArgs():
    parser = argparse.ArgumentParser(description=
                                     "Error rate estimations for context and"
                                     "cycle based errors")
    subparsers = parser.add_subparsers(dest='mode', help='sub-command help')
    parser_cycle = subparsers.add_parser('cycleError', help="blah")
    parser_cycle.add_argument("mdBam", help="bam generated by Samtools calmd"
                                            "command")
    parser_cycle.add_argument("-cycleheat", default=False, action='store_true')
    parser_cycle.add_argument("-family_size", default=None)
    parser_heatmap = subparsers.add_parser('heatmap', help='blah')
    parser_heatmap.add_argument("reads2", help="non-consensus fastq of read2")
    parser_heatmap.add_argument("bam", help="alignment of consensus reads")
    parser_heatmap.add_argument("kmerLen")
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


def errorFinder(read, readErr, readObs):
    seq = read.query_sequence
    start = read.qstart
    end = read.qend
    for base in xrange(start,end):
        readObs[base] += 1
        if seq[base] != '=' and seq[base] != 'N':
            readErr[base] += 1


def cycleError(args):
    rlen = pysam.AlignmentFile(args.mdBam).next().inferred_length
    #rlen = 146
    read1error = np.zeros(rlen,dtype=np.int)
    read1obs = np.zeros(rlen,dtype=np.int)
    read2error = np.zeros(rlen,dtype=np.int)
    read2obs = np.zeros(rlen,dtype=np.int)
    qc = 0
    fmc = 0
    rc = 0
    fam_range = args.family_size
    if fam_range:
        fam_range = [int(i) for i in fam_range.split(',')]
    for read in pysam.AlignmentFile(args.mdBam):
        if(read.is_secondary or read.is_supplementary or
           read.is_unmapped or read.is_qcfail):
            qc += 1
            continue
        if fam_range:
            if(read.opt('FM') < fam_range[0] or read.opt('FM') > fam_range[1]):
                fmc += 1
                continue
        if read.is_read1:
            rc += 1
            errorFinder(read, read1error, read1obs)
        if read.is_read2:
            rc += 1
            errorFinder(read, read2error, read2obs)
    sys.stdout.write("Family Size Range: %i-%i\n" % (fam_range))
    sys.stdout.write("Reads Analyzed: %i\n" % (rc))
    sys.stdout.write("Reads QC Filtered: %i\n" % (qc))
    if fam_range:
        sys.stdout.write("Reads Family Size Filtered: %i\n" % (fmc))
    read1prop = read1error/read1obs
    read2prop = read2error/read2obs
    read1mean = np.mean(read1prop)
    read2mean = np.mean(read2prop)
    sys.stdout.write("cycle\tread1\tread2\tread count\n")
    for index in xrange(rlen):
        sys.stdout.write("%i\t%f\t%f\t%i\n" % (index+1, read1prop[index],
                                             read2prop[index], rc))
    sys.stdout.write("%i\t%f\t%f\t%i\n" % (fam_range[1], read1mean, read2mean,
                                           rc))
    if args.cycleheat:
        cycleHeater(read1prop, read2prop, rlen)


def cycleHeater(read1prop, read2prop, rlen):
    fig, ax = plt.subplots()
    data = np.array([read1prop, read2prop])
    heatmap = ax.pcolor(data, cmap=plt.cm.Reds)
    colorb = plt.colorbar(heatmap)
    plt.show()

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
            raise ThisIsMadness("Did not find both reads in bam")
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
    args = getArgs()
    if args.mode not in ["heatmap", "cycleError"]:
        raise ThisIsMadness("invalid mode")
    if args.mode == "heater":
        heater(args)
    if args.mode == "cycleError":
        cycleError(args)


if __name__ == "__main__":
    main()
