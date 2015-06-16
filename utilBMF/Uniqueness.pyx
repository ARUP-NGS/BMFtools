"""
Contains a set of utilities for evaluating uniqueness of both
mappability and capture. Also might have applications in assembly-based
variant-calling.
"""
import cython
import pysam
import uuid
import logging
from subprocess import check_output, check_call
from utilBMF.HTSUtils import (GetUniqueItemsL, GetUniqueItemsD,
                              printlog as pl, ParseBed,
                              hamming_cousins, RevCmp)
from utilBMF.ErrorHandling import ThisIsMadness
from cytoolz import frequencies as cyfreq
from itertools import chain
from functools import partial
from operator import attrgetter
try:
    import re2 as re
except ImportError:
    import re
# Dictionary of edit distance flags to the reference (excludes clipping)
mmDict = {i: re.compile("NM:i:[0-%s]" % i) for i in xrange(20)}
oagseq = attrgetter("seq")
cfi = chain.from_iterable
hammingPt = partial(hamming_cousins, n=1)

cimport cython


cpdef cystr SequenceToFakeFq(cystr seq):
    return ("@" + seq + "\n" + seq +
            "\n+\n" + "G" * len(seq))


cpdef list GetKmersToCheck(cystr ref, int k=30, list bedline=[],
                           int padding=-1):
    """
    Gets a list of kmers which provide unique mappability
    to the region of interest.
    bedline should be the first 3 columns from a line in a bed file.
    """
    cdef int i, start, end
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


cpdef cystr FastqStrFromKmerList(list kmerList):
    """
    Creates a dummy fastq string from a list of kmers.
    """
    return "\n".join(map(SequenceToFakeFq, kmerList))


cdef class RefKmer(object):
    """
    Contains useful information regarding representative kmers selected
    from reference sequence.
    :param cystr seq
    :param cystr contig=None
    :param int pos=-1
    """

    def __init__(self, cystr seq, cystr contig=None,
                 int pos=-1):
        assert pos >= 0  # pos needs to be set
        self.contig = contig
        self.seq = seq
        self.len = len(seq)
        self.pos = pos

    @cython.returns(cystr)
    def __str__(self):
        return "%s|%s|%s" % (self.contig, self.pos, self.seq)


cdef class KmerFetcher(object):
    """Contains parameters for finding representative kmers.
    I want to permit the mixing of kmers of different sizes, so I have
    set a baseK, which is the K that we start using.

    minMQ defaults to 1, meaning a unique mapping for most aligners.
    Bowtie instead uses 255 to refer to unique alignments and (0,3) to mark
    multiple acceptable alignments.

    TODO: I'd like to have it know where dropout regions for these kmers
    are and have it increase the size of k only as necessary... Not sure
    how to automate that.

    :param cystr ref - path to reference fasta file.
    :param int padding - distance around region to pad for
           looking for kmers
    :param int mismatches - maximum permitted mismatches in alignment.
    Defaults to 0.
    :param int minMQ - minimum MQ for a kmer's alignment to be
    considered acceptable.
    Defaults to 1
    :param int k - length of kmers;
    * Note: minMQ doesn't make sense with bowtie/bowtie2's mapping quality
    assignments.
    Example use for creating a mapper with this:

    KF = KmerFetcher(ref="/home/daniel/human_g1k_v37.fasta", k=32)
    mapper = KF.BedToMap("~/cfDNA_targets.bed")

    Now you can map reads with this mapper!

    Alternatively, you can use KF.IBedToMap("~/cfDNA_targets.bed").
    Instead of assigning that hashmap to mapper as above, the hashmap
    is saved as self.FullMap.

    """
    def __init__(self, cystr ref=None, int padding=50,
                 int mismatches=0, int minMQ=1,
                 int k=30, int seed=-1, cystr aligner='mem'):
        self.ref = ref
        if(mismatches >= 0):
            self.mismatches = mismatches
        else:
            pl("KmerFetcher's mismatches set to below 0 - setting to 0, "
               "as that breaks downstream steps.")
            self.mismatches = 0
        self.minMQ = minMQ
        self.padding = padding
        self.k = k
        self.seed = seed
        self.aligner = aligner
        self.HashMap = {}
        self.FullMap = None

    def IBedToMap(self, cystr bedpath):
        """Fills the HashMap for each line in the bed file.
        Sets the FullMap hashmap to map all strings within hamming distance
        of mismatches limit.
        I is for in-place.
        """
        cdef list bedline
        cdef list bedlines = ParseBed(bedpath)
        for bedline in bedlines:
            self.FillMap(bedline)
        self.FullMap = self.BuildMapper(maplimit=self.mismatches)

    @cython.returns(dict)
    def BedToMap(self, cystr bedpath):
        """Fills the HashMap for each line in the bed file.
        Returns a hashmap which maps all strings within hamming distance
        of mismatches limit.
        """
        cdef list bedline
        cdef list bedlines = ParseBed(bedpath)
        for bedline in bedlines:
            self.FillMap(bedline)
        return self.BuildMapper(maplimit=self.mismatches)
    # Should we consider changing the number of mapping mismatches permitted
    # so as to be more stringent than with our mappability analysis?

    cdef setK(self, int newK):
        """Sets the value of K to the argument.
        """
        self.k = newK

    cdef int getK(self):
        return self.k

    cpdef cystr getFastqString(self, list bedline):
        return FastqStrFromKmerList(GetKmersToCheck(self.ref, k=self.k,
                                                    bedline=bedline,
                                                    padding=self.padding))

    cpdef cystr getOutputString(self, list bedline, cystr aligner="mem"):
        if(aligner == "mem"):
            return BwaFqToStr(self.getFastqString(bedline), ref=self.ref)
        elif(aligner=="bwt"):
            return BowtieFqToStr(self.getFastqString(bedline), ref=self.ref,
                                 seed=self.seed, mismatches=self.mismatches)
        else:
            raise ValueError("Sorry, only bwa mem or bowtie is supported currently.")

    cpdef FillMap(self, list bedline):
        self.HashMap[
            ":".join(map(str, bedline))] = self.GetUniqueKmers(bedline)

    cpdef list GetUniqueKmers(self, list bedline):
        return GetMQPassRefKmersMem(self.getOutputString(bedline, aligner=self.aligner),
                                    maxNM=self.mismatches, minMQ=self.minMQ)

    cpdef FMfrombed(self, cystr bedfile):
        cdef list bedline
        for bedline in ParseBed(bedfile):
            self.FillMap(bedline)

    def __getitem__(self, key):
        return self.HashMap[key]

    def iteritems(self, *args, **kwargs):
        return self.HashMap.iteritems(*args, **kwargs)

    def items(self, *args, **kwargs):
        return self.HashMap.items(*args, **kwargs)

    def itervalues(self, *args, **kwargs):
        return self.HashMap.itervalues(*args, **kwargs)

    def values(self, *args, **kwargs):
        return self.HashMap.values(*args, **kwargs)

    def iterkeys(self, *args, **kwargs):
        return self.HashMap.iterkeys(*args, **kwargs)

    def keys(self, *args, **kwargs):
        return self.HashMap.keys(*args, **kwargs)

    @cython.returns(dict)
    def BuildMapper(self, object cfi=cfi, int maplimit=1,
                    object oagseq=oagseq, object hammingPt=hammingPt,
                    object RevCmp=RevCmp):
        """
        Once all regions of interest have been build,
        create the mapping dictionary.
        :param cfi - localizes the chain.from_iterable call.
        :param hammingPt - localizes a partial hamming distance.
        """
        cdef cystr region, kmer
        cdef list kmerlist, mappingTups
        cdef dict mappingDict
        mappingTups = []
        for region, kmerlist in self.iteritems():
            mappingTups += [(kmer, region) for
                            kmer in list(cfi(map(hammingPt,
                                                 map(oagseq, kmerlist))))]
        return dict(mappingTups + [(RevCmp(kmer), region) for
                                   kmer, region in mappingTups])


@cython.returns(dict)
def GetRepresentativeKmerDict(*args, **kwargs):
    return cyfreq(GetRepKmersBwt(*args, **kwargs))


cpdef list GetRepKmersBwt(cystr ref, int k=30,
                          list bedline=[],
                          int padding=-1, int seedlen=-1,
                          int mismatches=-1, int minMQ=1,
                          bint useBowtie=False):
    cdef cystr fqStr, output
    fqStr = FastqStrFromKmerList(GetKmersToCheck(ref, k=k, bedline=bedline,
                                                 padding=padding))
    if(useBowtie):
        output = BowtieFqToStr(fqStr, ref=ref, seed=seedlen,
                               mismatches=mismatches)
        return GetUniqMQsBowtie(output, minMQ=minMQ)
    else:
        output = BwaFqToStr(fqStr, ref=ref, seed=seedlen)
    return GetUniqMQsBowtie(output, minMQ=minMQ)


cpdef cystr BowtieFqToStr(cystr fqStr, cystr ref=None,
                          int seed=-1, int mismatches=-1):
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
    outStr = check_output(cStr, shell=True)  # Capture output to string
    check_call(["rm", tmpFile])  # Delete the temporary file.
    print("Returning bowtieFqToStr output")
    return outStr


@cython.returns(cystr)
def BwaFqToStr(cystr fqStr, cystr ref=None,
               int seed=-1):
    """
    Returns the string output of a bwa mem call.
    """
    cdef cystr seedStr, cStr, outStr, tmpFile
    tmpFile = str(uuid.uuid4().get_hex()[0:8]) + ".hey.i.am.a.prefix.fq"
    tmpFileHandle = open(tmpFile, "w")
    tmpFileHandle.write(fqStr)
    tmpFileHandle.close()
    if(seed > 0):
        seedStr = "-k %s" % seed
    else:
        seedStr = ""
    cStr = "bwa mem -a %s %s %s" % (ref, seedStr, tmpFile)
    pl("Bowtie command string: %s" % cStr, level=logging.DEBUG)
    print("Bwa command string: %s" % cStr)
    outStr = check_output(cStr, shell=True)  # Capture output to string
    check_call(["rm", tmpFile])  # Delete the temporary file.
    print("Returning BwaFqToStr output")
    return outStr


@cython.returns(bint)
def PassesNM1(cystr rStr, int maxNM=2, dict mmDict=mmDict):
    """
    Checks a SAM line to see if its edit distance is below or equal
    to the maximum.
    """
    cdef int i
    cdef list strList = ["NM:i:%s" % i for i in xrange(maxNM + 1)]
    for item in strList:
        if item in rStr:
            return True
    return False


@cython.returns(bint)
def PassesNM(cystr rStr, int maxNM=2, dict mmDict=mmDict):
    """
    Checks a SAM line to see if its edit distance is below or equal
    to the maximum.
    """
    return mmDict[maxNM].match(rStr) is not None


@cython.returns(list)
def GetMQPassRefKmersMem(cystr bwaStr, int maxNM=2, int minMQ=1):
    """
    Takes a string output from bowtie (SAM format) and gets the names of the reads
    with MQ >= minMQ that are unique alignments. Returns a RefKmer object
    built from passing bowtie output with unique mappings (no XA:Z or XT:A:R flags or
    non-zero mapping qualities).
    """
    cdef list lines, i
    cdef cystr f
    cdef tuple nameCount
    return [RefKmer(i[0], contig=i[2],
                    pos=int(i[3])) for i in [f.strip().split("\t") for
                                             f in bwaStr.split("\n") if
                                             "XA:Z:" not in f and  # Supp aln
                                             "XT:A:R" not in f and  # Supp aln
                                             f != "" and  # Empty line
                                             f[0] != "@" and  # Header
                                             PassesNM(f, maxNM=maxNM)]  # NM
            if i[4] != "0" and i[4] >= minMQ]


@cython.returns(list)
def GetMQPassReadsMem(cystr bwaStr):
    """
    Takes a string output from bowtie and gets the names of the reads
    with MQ >= minMQ. Defaults to 1 (for a unique alignment)
    """
    cdef list i
    cdef cystr f
    return [i[0] for i in [f.strip().split("\t") for
                           f in bwaStr.split("\n") if f != "" and
                           f[0] != "@" and "XA:Z:" not in f]
            if i[4] != "0"]


@cython.returns(list)
def GetMQPassRefKmersBwt1(cystr bwtStr):
    """
    Takes a string output from bowtie and gets the names of the reads
    with MQ >= minMQ. Defaults to 1 (for a unique alignment)
    """
    cdef list lines, i
    cdef cystr f
    cdef tuple nameCount
    return [RefKmer(i[0], contig=i[2], pos=int(i[3])) for
            i in [f.strip().split("\t") for
                  f in bwtStr.split("\n") if
                  f != "" and f[0] != "@"]
            if i[4] == "255"]


@cython.returns(list)
def GetMQPassReadsBwt1(cystr bwtStr):
    """
    Takes a string output from bowtie and gets the names of the reads
    with MQ >= minMQ. Defaults to 1 (for a unique alignment)
    """
    cdef list lines, i
    cdef cystr f
    cdef tuple nameCount
    return [i[0] for i in [f.strip().split("\t") for
                           f in bwtStr.split("\n") if
                           f != "" and f[0] != "@"]
            if i[4] == "255"]


@cython.returns(list)
def GetUniqMQsBowtie(cystr bwtStr, int minMQ=1):
    """
    Takes a string output from bowtie and gets the names of the reads
    with MQ >= minMQ. Defaults to 1 (for a unique alignment)
    """
    return GetUniqueItemsL(GetMQPassRefKmersBwt1(bwtStr))
