from __future__ import division
import numpy as np
from utilBMF.HTSUtils import pFastqProxy, pFastqFile, getBS
from itertools import groupby
from utilBMF.ErrorHandling import ThisIsMadness

consInFq = pFastqFile("/home/brett/Projects/BMFTools_Devel/lamda_data/kmer_t"
                        "est/lamda_test_R1.rescued.shaded.BS.cons.fastq")
readsInFq = pFastqFile("/home/brett/Projects/BMFTools_Devel/lamda_data/kmer_t"
                        "est/lamda_test_R1.rescued.shaded.BS.fastq")

for bc4fq, fqRecGen, in groupby(readsInFq, key=getBS):
    consRead = consInFq.next()
    consBS = getBS(consRead)
    if(consBS != bc4fq):
        raise ThisIsMadness("Barcode mismatch, are both read files sorted?")
    pFqPrxList = list(fqRecGen)
    for rec in pFqPrxList:
        if(rec.sequence != consRead.sequence):
            print rec.sequence, getBS(rec)
            print consRead.sequence, getBS(consRead)
            print
