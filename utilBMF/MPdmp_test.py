from MawCluster.BCFastq import FastFisherFlattening, MakeTagComment
from utilBMF.HTSUtils import pFastqProxy, pFastqFile, permuteNucleotides

import numpy as np
import multiprocessing as mp

import argparse
<<<<<<< HEAD
import pysam
=======


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
        """
        1. get all records with same barcode into dict
        2. Mark all records with BS sequence
        3. Sort records by BS sequence (groupby?)
        4. Consolidate records with the same sequence
        5. Print records to a temp file with the prefix sequence
        6. Delete fqRecords dictionary
        """
        pass

    def _finish(self):
        self.jobComplete=True
        self.fqRecords={}

>>>>>>> bbd9a28459ffcf98639bda160e42fdf9e99a5266

def shadeRead(read, BC, head=0):
    if("N" in BC):
        read.comment += "|FP=0|BS=%s" % (read.sequence[1:head + 1] + BC)
    else:
        read.comment += "|FP=1|BS=%s" % (read.sequence[1:head + 1] + BC)
    return read

def worker(Fastq, IndexFq, bcLen, prefix, head=0):
    """
    class based implementation above probably won't work.  switching
    to a function-based one, using a worker function.  I should probably
    have one job that queues workers then a worker function that actually
    executes, going to need to compare pool vs queue.
    """
    fastq = pFastqFile(Fastq)
    indexFq = pFastqFile(IndexFq)
    bcHash = {}
    ifn = indexFq.next
    lenPrefix = len(prefix)
    print("now starting prefix %s" % (prefix))
    for read in fastq:
        BC = ifn().sequence
        if prefix == BC[:lenPrefix]:
            try:
                bcHash[BC].append(shadeRead(read, BC))
            except KeyError:
                bcHash[BC] = [shadeRead(read, BC)]
    fqname = Fastq.split('.')[0]
    output = open(fqname+"."+prefix+".fastq",'w')
    for barcode in bcHash.iterkeys():
        output.write(FastFisherFlattening(bcHash[barcode], barcode))
    output.close()
    numBC = len(bcHash)  # returns 0 if uninitialized, length otherwise.
    del bcHash
    print("prefix %s complete" % (prefix))
    return numBC


def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("Fastq1", help="Fastq1", type=str)
    parser.add_argument("Fastq2", help="Fastq2", type=str)
    parser.add_argument("IndexFq", help="Barcode Fastq", type=str)
    parser.add_argument("ncpus", help="number of cpus", type=int)
    parser.add_argument("lenPrefix", help="hacky, temporary length of prefix"
                        " by which we will split the jobs",
                        type=int)
    return parser.parse_args()

if __name__ == '__main__':
    args = getArgs()
    pool = mp.Pool(processes=args.ncpus)
<<<<<<< HEAD
    bcLen = len(pFastqFile(args.IndexFq).next().sequence)
    allPrefixes = permuteNucleotides(4**args.lenPrefix,kmerLen=args.lenPrefix)
    results = [pool.apply_async(worker, args=(args.Fastq1, args.IndexFq, bcLen,
                p)) for p in allPrefixes]
=======
    allPrefixes = permuteNucleotides(4**args.lenPrefix, kmerLen=args.lenPrefix)
    results = [pool.apply_async(worker, args=(args.Fastq1, args.IndexFq, p, ))
               for p in allPrefixes]
>>>>>>> bbd9a28459ffcf98639bda160e42fdf9e99a5266
    things = [p.get() for p in results]
    print things
