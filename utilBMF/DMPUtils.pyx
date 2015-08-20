import logging
import gzip
import os
import subprocess
import time

from MawCluster.BCFastq import FastFisherFlattening, CutadaptPaired
from utilBMF.HTSUtils import (pFastqProxy, pFastqFile, permuteNucleotides,
                                printlog as pl, getBS)
from utilBMF.ErrorHandling import ThisIsMadness as Tim
from utilBMF.ErrorHandling import UnsetRequiredParameter, ImproperArgumentError

import multiprocessing as mp

def singleDMPWorker(Fastq, IndexFq, bcLen, Prefix, hpLimit, Head=0):
    """
    Not implemented!
    """
    raise Tim("Not yet Implemented!  Sorry!")
    pass


def pairedDMPWorker(Fastq1, Fastq2, IndexFq, bcLen, Prefix, hpLimit, Head=0,
                    cutAdapt=False, p3Seq="default", p5Seq="default"):
    """
    dmpWorker function, pulls a bunch of reads into memory based on their
    BC prefix marks then consolidates them
    """
    cdef int lenPrefix
    cdef pFastqFile_t fq1, fq2, indexFq
    cdef pFastqProxy_t read1, read2
    cdef cystr BC, bin
    cdef dict bcHash1, bcHash2
    fq1 = pFastqFile(Fastq1)
    fq2 = pFastqFile(Fastq2)
    indexFq = pFastqFile(IndexFq)
    ifn = indexFq.next
    ifq2n = fq2.next
    bcHash1 = {}
    bcHash2 = {}
    lenPrefix = len(Prefix)
    pl("now starting dmp on prefix %s" % (Prefix))
    for read1 in fq1:
        BC = ifn().sequence
        read2 = ifq2n()
        if Head:
            BC = (read1.sequence[1:Head + 1] + BC +
                       read2.sequence[1:Head + 1])
        bin = BC[:lenPrefix]
        if Prefix == bin:
            read1.comment = cMakeTagComment(BC, read1, hpLimit)
            read2.comment = cMakeTagComment(BC, read2, hpLimit)
            try:
                bcHash1[BC].append(read1)
                bcHash2[BC].append(read2)
            except KeyError:
                bcHash1[BC]=[read1]
                bcHash2[BC]=[read2]
    fq1name = Fastq1.split(".")[0]+"."+Prefix+".fastq"
    fq2name = Fastq2.split(".")[0]+"."+Prefix+".fastq"
    output1 = open(fq1name,'w')
    output2 = open(fq2name,'w')
    with open(fq1name, 'w') as output1, open(fq2name, 'w') as output2:
        for barcode in bcHash1.keys():
            output1.write(FastFisherFlattening(bcHash1[barcode], barcode))
            output2.write(FastFisherFlattening(bcHash2[barcode], barcode))
    if cutAdapt:
        cafq1, cafq2 = CutadaptPaired(fq1name, fq2name, p3Seq, p5Seq)
        pl("runnning cutadapt on temporary fastqs with prefix %s" %(Prefix))
        subprocess.check_call("rm", fq1name)
        subprocess.check_call("rm", fq2name)
    pl("completed dmp on prefix %s" % (Prefix))
    return fq1name, fq2name


def get_uncompressed_size(file):
    """
    Used to get uncompressed Fastq file sizes
    possibly ineffiecent, not currently Used
    """
    fileobj = open(file, 'r')
    fileobj.seek(-8, 2)
    crc32 = gzip.read32(fileobj)
    isize = gzip.read32(fileobj)  # may exceed 2GB
    fileobj.close()
    return isize


def calcMaxRam(maxRam, ncpus, fastq):
    """
    Someday I will use this to estimate the proper prefix length
    for splitting a fastq based on desired maxRam and nCPUs.
    Today is not that day.
    """
    cdef long fqSize, mRam
    cdef int lprefix
    mRam = maxRam*1000000000
    fqSize = get_uncompressed_size(fastq)*2
    for x in range(3,10):
        topRam = (fqSize/(4**x))*ncpus
        if topRam < mRam:
            lprefix = x
            break
    return lprefix


def multiProcessingDemulitplex(inFqs, bcLen, indexFq="default", int head=-1,
                                int ncpus=1, int num_nucs=3, hpLimit=1):
    """
    Args I need:
        1. In Fastq, 1 or 2 fastqs in a list, call worker based on number
        2. Index fastq containing barcodes
        3. amount of sequence to salt into barcode (head) default = 0
        4. homopolymer proportion (default 80%) to fail reads
        5. Runnign cut adapt? True/False
        6. if above true, p3Seq, p5Seq, and overlap len
        7. number of CPUS to run across
        8.
    Similar to DispatchParallelDMP, but redesigned for multiprocessing
    of the DMP process, however, here we do not call to cutadapt,
    that will be handled downstream.
    """
    if num_nucs != 3:
        pl("using custom number of nucleotides for splitting temporary files"
        ", recommended value is 3, used here %s" % (num_nucs))
    if indexFq == "default":
        raise UnsetRequiredParameter("Index fastq with barcodes is required")
    if head < 0:
        raise UnsetRequiredParameter("Improper or unset head value defaulting"
                                     "to 0")
        head = 0
    pl("running multiprocessed DMP using %s CPUs" % (ncpus))
    allPrefixe = permuteNucleotides(4**num_nucs, kmerLen=num_nucs)
    pass
