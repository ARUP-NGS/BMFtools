"""
Contains a set of utilities for evaluating uniqueness of both
mappability and capture. Also might have applications in assembly-based
variant-calling.
"""
import cython
import pysam
import uuid
import logging
from subprocess import check_output
from utilBMF.HTSUtils import GetUniqueItemsL, GetUniqueItemsD, ThisIsMadness, printlog as pl
from cytoolz import frequencies as cyfreq
cimport cython

@cython.returns(cython.str)
def SequenceToFakeFq(cython.str seq):
    return ("@" + seq + "\n" + seq +
            "\n+\n" + "G" * len(seq))


@cython.returns(list)
def GetKmersToCheck(cython.str ref, cython.int k=30, list bedline=[],
                    cython.int padding=-1):
    """
    Gets a list of kmers which provide unique mappability
    to the region of interest.
    bedline should be the first 3 columns from a line in a bed file.
    """
    cdef cython.int i, start, end
    cdef list kmerList, candidateKmers
    if(padding < 0):
        pl("Padding not set - defaults to kmer size.")
        padding = k
    kmerList = []
    refHandle = pysam.FastaFile(ref)
    contig, start = bedline[0], bedline[1] - padding
    end = bedline[2] + padding
    regionStr = refHandle.fetch(contig, start, end)
    return [regionStr[i:i + k] for i in xrange(end - start - k)]


@cython.returns(cython.str)
def FastqStrFromKmerList(list kmerList):
    """
    Creates a dummy fastq string from a list of kmers.
    """
    return "\n".join(map(SequenceToFakeFq, kmerList))


cdef class KmerFetcher(object):
    """
    Contains parameters for finding representative kmers.
    I want to permit the mixing of kmers of different sizes, so I have
    set a baseK, which is the K that we start using.
    """
    def __init__(self, cython.str ref=None, cython.int padding=-1,
                 cython.int mismatches=-1, cython.int minMQ=1,
                 cython.int k=30):
        self.ref = ref
        self.mismatches = mismatches
        self.minMQ = minMQ
        self.padding = padding
        self.k = k
        self.HashMap = {}

    cdef setK(self, cython.int newK):
        self.k = newK

    cdef cython.int getK(self):
        return self.k

    cdef cython.str GetBowtieOutput(self, cython.str fqStr):
        return BowtieFqToStr(fqStr, ref=self.ref,
                             seed=self.k, mismatches=self.mismatches)

    cdef FillMap(self, list bedline):
        self.HashMap[":".join(map(str, bedline))] = self.GetUniqueKmers(bedline)

    cdef list GetUniqueKmers(self, list bedline):
        return GetUniqueItemsL(self.GetKmerList(bedline))

    cdef list GetKmerList(self, list bedline):
        return GetRepresentativeKmerList(self.ref, k=self.k, bedline=bedline,
                                         mismatches=self.mismatches,
                                         minMQ=self.minMQ,
                                         padding=self.padding, seedlen=self.k)

@cython.returns(dict)
def GetRepresentativeKmerDict(*args, **kwargs):
    return cyfreq(GetRepresentativeKmerList(*args, **kwargs))


@cython.returns(list)
def GetRepresentativeKmerList(cython.str ref, cython.int k=30,
                              list bedline=[],
                              cython.int padding=-1, cython.int seedlen=-1,
                              cython.int mismatches=-1, cython.int minMQ=1):
    cdef cython.str fqStr, btOutStr
    fqStr = FastqStrFromKmerList(GetKmersToCheck(ref, k=k, bedline=bedline,
                                                 padding=padding))
    btOutStr = BowtieFqToStr(fqStr, ref=ref, seed=seedlen, mismatches=mismatches)
    return GetUniqMQs(btOutStr, minMQ=minMQ)


@cython.returns(cython.str)
def BowtieFqToStr(cython.str fqStr, cython.str ref=None,
                  cython.int seed=-1, cython.int mismatches=-1):
    """
    Returns the string output of a bowtie2 call.
    With bowtie, you can specify precisely the number of permitted mismatches
    in the string, but with bowtie2, you don't need to write any temporary
    files. I'll play around with the alignment settings as I get along in
    the getkmers tool.
    """
    if(seed < 0):
        raise ThisIsMadness("seed length must be set for Bowtie2FqToStr.")
    if(mismatches > 2 or mismatches < 0):
        raise ThisIsMadness("0 <= mismatches <= 2!")
    tmpFile = str(uuid.uuid4().get_hex()[0:8]) + ".hey.i.am.a.prefix.fq"
    tmpFileHandle = open(tmpFile, "w")
    tmpFileHandle.write(fqStr)
    tmpFileHandle.close()
    cStr = "bowtie --mm --all -n %s -l %s %s -S %s" % (mismatches, seed,
                                                       ref, tmpFile)
    pl("Bowtie command string: %s" % cStr, level=logging.DEBUG)
    print("Bowtie command string: %s" % cStr)
    outStr = check_output(cStr, shell=True) # Capture output to string
    # check_call(["rm", tmpFile])  # Delete the temporary file.
    print("Returning bowtieFqToStr output")
    return outStr


@cython.returns(list)
def GetMQPassReadsBowtie1(cython.str bwtStr):
    """
    Takes a string output from bowtie and gets the names of the reads
    with MQ >= minMQ. Defaults to 1 (for a unique alignment)
    """
    cdef list lines, i
    cdef cython.str f
    cdef tuple nameCount
    return [i[0] for i in [f.strip().split("\t") for
                           f in bwtStr.split("\n") if
                           f != "" and f[0] != "@"]
            if i[4] == "255"]


@cython.returns(list)
def GetUniqMQs(cython.str bwtStr, cython.int minMQ=1):
    """
    Takes a string output from bowtie and gets the names of the reads
    with MQ >= minMQ. Defaults to 1 (for a unique alignment)
    """
    return GetUniqueItemsL(GetMQPassReadsBowtie1(bwtStr))