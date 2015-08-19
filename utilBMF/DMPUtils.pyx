
import logging
import gzip
import time

from MawCluster.BCFastq import FastFisherFlattening, MakeTagComment
from utilBMF.HTSUtils import (pFastqProxy, pFastqFile, permuteNucleotides,
                                printlog as pl, getBS)
from utilBMF.ErrorHandling import ThisIsMadness as Tim
from utilBMF.ErrorHandling import UnsetRequiredParameter, ImproperArgumentError

import multiprocessing as mp

def shadeRead(read, BC):
    pass

def dmpWorker(Fastq, IndexFq, bcLen, Prefix, Head=0):
    """
    dmpWorker function, pulls a bunch of reads into memory based on their
    BC prefix marks then consolidates them
    """
    cdef int lenPrefix
    # cdef pFastqFile_t fastq, indexFq

    fastq = pFastqFile(Fastq)
    indexFq = pFastqFile(IndexFq)
    ifn = indexFq.next
    bcHash = {}
    lenPrefix = len(Prefix)

    pl("now starting dmp on prefix %s" % (Prefix))

    for read in fastq:
        BC = ifn().sequence
        if Prefix == BC[:lenPrefix]:
            try:
                bcHash[BC].append(shadeRead(read, BC))
            except KeyError:
                bcHash[BC]=[shadeRead(read, BC)]
    fqname = Fastq.split(".")[0]+"."+Prefix+".fastq"
    output = open(fqname,'w')
    for barcode in bcHash.keys():
        output.write(FastFisherFlattening(bcHash[barcode], barcode))
    output.close()
    pl("completed dmp on prefix %s" % (Prefix))
    return fqname

def multiProcessingDemulitplex(inFqs, indexFq="default", int head=-1,
                                int ncpus=1, int num_nucs=3):
    """
    Similar to DispatchParallelDMP, but redesigned for multiprocessing
    of the DMP process, however, here we do not call to cutadapt,
    that will be handled downstream.
    """
    if num_nucs != 3:
        pl("using custom number of nucleotides for smaller, more numerous"
        "temporary files, recommended value is 3, used here %s" % (num_nucs))
    if indexFq == "default":
        raise UnsetRequiredParameter("Index fastq with barcodes is required")
    if head < 0:
        raise UnsetRequiredParameter("Improper or unset head value defaulting"
                                     "to 0")
        head = 0
    pass
