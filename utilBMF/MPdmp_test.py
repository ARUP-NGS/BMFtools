from utilBMF.HTSUtils import pFastqProxy, pFastqFile, permuteNucleotides

import numpy as np
import multiprocessing as mp
import pysam

import argparse

class shepard(object):
    """
    pushes stuff around
    """

    def __init__(self, fastq, indexFq, bsPrefix):
        self.prefix=bsPrefix
        self.prefixLen=len(bsPrefix)
        self.Fq=fastq
        self.iFq=indexFq
        self.jobComplete=False
        self.fqRecords={}

    def _run(self):
        print "running prefix: ", self.prefix
        for record in self.Fq:
            BSrec=self.iFq.next()
            BS=BSrec.seq
            if BS[self.prefixLen] == self.prefix:
                self.fqRecords[prefix].append(record)


def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("Fastq1", help="Fastq1", type=str)
    parser.add_argument("Fastq2", help="Fastq2", type=str)
    parser.add_argument("IndexFq", help="Barcode Fastq", type=str)
    parser.add_argument("ncpus", help="number of cpus", type=int)
    parser.add_argument("lenPrefix", help="hacky, temporary length of prefix"
                        " by which we will split the jobs, complicated.",
                        type=int)
    return parser.parse_args()

if __name__ == '__main__':
    args = getArgs()
    pool = mp.Pool(processes=args.ncpus)
    allPrefixes = permuteNucleotides(4**args.lenPrefix,kmerLen=args.lenPrefix)
    fastq1 = pFastqFile(args.Fastq1)
    fastq2 = pFastqFile(args.Fastq2)
    indexFq = pFastqFile(args.IndexFq)
    flock=[]
    for prefix in allPrefixes:
        flock.append(shepard(fastq1, indexFq, prefix))
        flock.append(shepard(fastq2, indexFq, prefix))
