# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True
import operator
import pysam
import cython
from cytoolz import map as cmap
from utilBMF.HTSUtils import (FractionAligned, FractionSoftClipped,
                              printlog as pl, CoorSortAndIndexBam,
                              GetKmersToCheck, FastqStrFromKmerList,
                              BowtieFqToStr, GetMQPassReads,
                              GetInsertionFromAlignedSegment,
                              GetDeletionFromAlignedSegment,
                              FacePalm, ParseBed,
                              shen, ssStringFromRead, ccopy)
from .SNVUtils import (HeaderInfoLine, HeaderFormatLine,
                       HeaderContigLine, HeaderCommandLine,
                       HeaderReferenceLine, HeaderFileFormatLine,
                       GetContigHeaderLines, HeaderFilterLine)
from utilBMF.ErrorHandling import ThisIsMadness
import os.path
import uuid
from cytoolz.itertoolz import frequencies as cyfreq
import sys

cimport pysam.calignmentfile
cimport cython
cimport numpy as np
cimport utilBMF.HTSUtils
ctypedef IndelQuiver IndelQuiver_t
ctypedef AbstractIndelContainer AbstractIndelContainer_t
ctypedef Insertion Insertion_t
ctypedef Deletion Deletion_t

"""
Contains utilities for working with indels for HTS data.
"""


def FilterByIndelRelevance(inBAM, indelOutputBAM="default",
                           otherOutputBAM="default",
                           minFamSize=2):
    """
    Writes reads potentially relevant to an indel to indelOutputBAM and
    other reads to the otherOutputBAM.
    idRel stands for indel relevant.
    idIrl stands for indel irrelevant.
    Input BAM must be name sorted: coordinate sorted is not supported!
    """
    cdef cython.long idrel
    cdef cython.long idirl
    cdef pysam.calignmentfile.AlignedSegment read1, read2, entry
    if(indelOutputBAM == "default"):
        indelOutputBAM = ".".join(inBAM.split(".")[0:-1] + ["indRel", "bam"])
    if(otherOutputBAM == "default"):
        otherOutputBAM = ".".join(inBAM.split(".")[0:-1] + ["indIrl", "bam"])
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    indelHandle = pysam.AlignmentFile(indelOutputBAM, "wb", template=inHandle)
    otherHandle = pysam.AlignmentFile(otherOutputBAM, "wb", template=inHandle)
    ohw = otherHandle.write
    ihw = indelHandle.write
    idrel = 0
    idirl = 0
    for entry in inHandle:
        if entry.is_read1:
            read1 = entry
            continue
        else:
            read2 = entry
        assert read1.query_name == read2.query_name
        if(IsIndelRelevant(read1, minFam=minFamSize) or
           IsIndelRelevant(read2, minFam=minFamSize)):
            ihw(read1)
            ihw(read2)
            idrel += 1
        else:
            ohw(read1)
            ohw(read2)
            idirl += 1
    inHandle.close()
    otherHandle.close()
    indelHandle.close()
    pl("Finished filtering by indel relevance. Relevant pairs: %s" % idrel +
       ". Irrelevant pairs: %s" % idirl)
    return indelOutputBAM, otherOutputBAM


def IsIndelRelevant(
        pysam.calignmentfile.AlignedSegment read, cython.long minFam=2,
        cython.float minSF=0.5, cython.bint keepUnmapped=False):
    """
    True if considered relevant to indels.
    False otherwise.
    """
    if(read.opt("FM") < minFam):
        return False
    if(read.cigarstring is None):
        # This read is simply unmapped. Let's give it a chance!
        if(keepUnmapped):
            return True
        return False
    if("I" in read.cigarstring or "D" in read.cigarstring):
        return True
    if(FractionSoftClipped(read) >= minSF):
        return True
    return False


def GetSCFractionArray(inBAM):
    cdef pysam.calignmentfile.AlignedSegment i
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    return [FractionSoftClipped(i) for i in inHandle]


def GetFreebayesCallStr(inBAM, ref="default", bed="default",
                        outVCF="default", ploidy=-1,
                        minMQ=1, minBQ=3, haplotypeLength=80,
                        rdf=1.0, bestNAlleles=100):
    """
    Used to call freebayes for indel calling.
    """
    if(os.path.isfile(inBAM + ".bai") is False):
        pl("BAM must be coordinate-sorted and "
           "indexed for freebayes to call variants. Doing so!")
        inBAM = CoorSortAndIndexBam(inBAM)
    if(outVCF == "default"):
        outVCF = ".".join(inBAM.split(".")[:-1]) + ".fb.vcf"
    if(ref == "default"):
        raise ThisIsMadness("Reference required for freebayes call.")
    if(bed == "default"):
        raise ThisIsMadness("Bed file require for freebayes call currently. "
                            "If there is sufficient demand for an all-regions"
                            " call, I can add that option.")
    if(ploidy < 0):
        pl("Ploidy less than 0. Defaulting to report best N alleles.")
        cStr = ("freebayes -v %s -t %s -f %s " % (outVCF, bed, ref) +
                "-q %s -m %s --haplotype-length" % (minBQ, minMQ) +
                " %s -D %s " % (haplotypeLength, rdf) +
                "--min-alternate-fraction 0 --pooled-continuous" +
                " %s" % inBAM)
    else:
        pl("Ploidy set: %s" % ploidy)
        cStr = ("freebayes -v %s -t %s -f %s " % (outVCF, bed, ref) +
                "-q %s -m %s --haplotype-length" % (minBQ, minMQ) +
                " %s -D %s " % (haplotypeLength, rdf) +
                "--min-alternate-fraction 0 --pooled-continuous" +
                " -p %s --use-best-n-alleles %s" % (ploidy, bestNAlleles) +
                " %s" % inBAM)
    return cStr


@cython.returns(cython.str)
def GetFBOutVCFFromStr(cython.str cStr):
    """
    Gets out vcf from freebayes call. Used for parallelization.
    """
    return cStr.split(" ")[2]


@cython.returns(list)
def GetUniquelyMappableKmers(cython.str ref, cython.long k=30,
                             list bedline=[], cython.long minMQ=1,
                             cython.long padding=-1, cython.long mismatches=2):
    """
    Uses a set of HTSUtils methods to find kmers from a region
    which are uniquely mappable. This makes it possible to do alignment-free
    variant-calling. (Well, except for the bwasw step).
    If no outfile is specified, defaults to stdout
    """
    cdef list kmerList, PassingReadNames
    cdef cython.str fqStr, bowtieStr
    pl("Getting potential kmers")
    kmerList = GetKmersToCheck(ref, k=k, bedline=bedline, padding=padding)
    pl("Making dummy fastq records for each kmer")
    fqStr = FastqStrFromKmerList(kmerList)
    pl("Aligning these kmers to the genome to test for unique mappability"
       " with a given number of mismatches %s and minMQ %s." % (mismatches,
                                                                minMQ))
    bowtieStr = BowtieFqToStr(fqStr, ref=ref, mismatches=mismatches, seed=k)
    PassingReadNames = GetMQPassReads(bowtieStr, minMQ=minMQ)
    return PassingReadNames


@cython.returns(IndelQuiver_t)
def FillIndelQuiverRegion(inBAM, cython.long minPairs=2,
                          cython.float minShen=0.2,
                          cython.long window=16, cython.str ref=None,
                          list bedRegion=[], cython.long minMQ=1,
                          cython.long minFM=1, cython.long minNumSS=1):
    """
    To fill an IndelQuiver, the bam must have the SV tags set, and positive
    calls require a DSI or DSD tag.

    You do need to be somewhat careful about adding Insertion or Deletion
    objects to an IndelQuiver object, as those objects could have been made
    with reads not satisfying minMQ and minFM.
    """
    cdef IndelQuiver_t Quiver
    cdef pysam.calignmentfile.AlignedSegment rec
    cdef pysam.calignmentfile.AlignmentFile inHandle
    cdef cython.long end
    Quiver = IndelQuiver(ref=ref, window=window, minFM=minFM, minMQ=minMQ,
                         bam=inBAM, minShen=minShen, minPairs=minPairs,
                         minNumSS=minNumSS)
    refHandle = pysam.FastaFile(ref)
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    inFetch = inHandle.fetch(bedRegion[0], bedRegion[1])
    end = bedRegion[2]
    for rec in inFetch:
        if(rec.mapping_quality < minMQ):
            continue
            """
            Since we're looking for DS[ID], if one read has a mapq
            sufficient for inclusion, we get the benefit of seeing that indel
            even if one of the reads failed the MQ threshold.
            """
        if(rec.opt("FM") < minFM):
            continue
            """
            Only working with read families with FM >= minFM
            cleans up our results.
            """
        try:
            svtags = rec.opt("SV")
        except KeyError:
            FacePalm("SV tags must be marked to call tagged indels!")
        if("DSI" in svtags):
            Quiver.addIndel(GetInsertionFromAlignedSegment(rec,
                                                           handle=refHandle))
        if("DSD" in svtags):
            Quiver.addIndel(GetDeletionFromAlignedSegment(rec,
                                                          handle=refHandle))
        if(rec.reference_pos >= end):
            pl("Finished bed region %s." % bedRegion)
            break
    return Quiver


@cython.returns(IndelQuiver_t)
def FillEntireQuiver(inBAM, cython.long minPairs=2, cython.float minShen=0.2,
                     cython.long window=16, cython.str ref=None,
                     cython.str bed=None, cython.long minFM=1,
                     cython.long minMQ=1, cython.long minNumSS=1):
    """
    Simply fills a quiver for each region, merges the quivers, and returns
    the full quiver.
    """
    cdef IndelQuiver_t Quiver, RegionalQuiver
    cdef list bedlines, line
    Quiver = IndelQuiver(ref=ref, window=window, minFM=minFM, minMQ=minMQ,
                         bam=inBAM, minShen=minShen, minPairs=minPairs,
                         minNumSS=minNumSS)
    bedlines = ParseBed(bed)
    for line in bedlines:
        RegionalQuiver = FillIndelQuiverRegion(inBAM, window=window,
                                               minFM=minFM, minMQ=minMQ,
                                               minPairs=minPairs,
                                               minShen=minShen, ref=ref,
                                               bedRegion=line,
                                               minNumSS=minNumSS)
        Quiver.mergeQuiver(RegionalQuiver)
    return Quiver


cdef class AbstractIndelContainer(object):
    """
    Base class for insertion and deletion container objects.

    Type can be -1, 0, or 1. (Deletion, deletion and insertion, and
    just insertion)
    Start and end refer to different things for insertions and deletions.
    For a deletion, start is the first missing base and end is the last missing
    reference base position.
    seq should be None for a deletion
    """

    def __init__(self, cython.str contig, cython.long start=-666,
                 cython.long end=-1, cython.long type=-137,
                 cython.str seq=None):
        self.contig = contig
        self.start = start
        self.end = end
        self.type = type
        self.seq = seq
        self.readnames = []
        self.uniqStr = None
        self.StartStops = []

    def __str__(self):
        raise ThisIsMadness("Abstract method must be inherited. Sorry, cdef w"
                            "on't let me actually make this an abstract "
                            "class.")

    @cython.returns(cython.long)
    def __len__(self):
        """
        Returns the number of reads supporting it which have been queried
        against this object since its creation.
        """
        return len(self.readnames)

    @cython.returns(cython.bint)
    def compare(self, indelObj):
        """
        Used for comparing two insertion or deletion objects.
        This just makes it a little cleaner and more flexible.
        That way, if they have the same unique identifier string
        but a different number of readnames.
        """
        assert isinstance(indelObj, AbstractIndelContainer)
        return self.uniqStr == indelObj.uniqStr

    def add(self, indelObj, inplace=True):
        """
        If the uniqStr attributes are identical, merge the families.
        I imagine that this would be useful when comparing lots of sets
        of indels across different samples or doing multiple passes.
        (e.g., duplex-supported indels and otherwise)
        """
        assert isinstance(indelObj, AbstractIndelContainer)
        if self.compare(indelObj):
            if(inplace):
                self.readnames += indelObj.readnames
                return self
            else:
                newIndelObj = ccopy(indelObj)
                newIndelObj.readnames += self.readnames
                return newIndelObj

    @cython.returns(dict)
    def GetNameCounter(self):
        return cyfreq(self.readnames)

    @cython.returns(cython.str)
    def __getitem__(self, cython.long index):
        return self.readnames[index]

    @cython.returns(list)
    def sort(self):
        self.readnames = sorted(self.readnames)

    @cython.returns(cython.long)
    def getNumSS(self):
        return(len(set(self.StartStops)))

    def register(self, pysam.calignmentfile.AlignedSegment read):
        self.readnames.append(read.query_name)
        self.StartStops.append(ssStringFromRead(read))

    def merge(self, AbstractIndelContainer_t AIC):
        try:
            assert AIC.uniqStr == self.uniqStr
        except AssertionError:
            raise ThisIsMadness("To merge two IndelContainer objects, "
                                "their unique string description must match.")
        self.readnames += AIC.readnames
        self.StartStops += AIC.StartStops


cdef class Insertion(AbstractIndelContainer):
    """
    Expansion of the AbstractIndelContainer object.
    "Start" is actually the preceding base to the insertion (counting up),
    "End" is the following. end should always be greater than start.
    The important thing is to be able to hold read names and have a unique
    string representing each indel so that we can make calls.

    If no handle is provided, shen (Shannon Entropy) is set to -1.
    """

    def __init__(self, pysam.calignmentfile.AlignedSegment read,
                 cython.str contig, cython.long start=-1,
                 cython.str seq=None, pysam.cfaidx.FastaFile handle=None,
                 cython.long window=20):
        if(start < 0):
            raise ThisIsMadness("start required for InsertionContainer.")
        self.contig = contig
        self.start = start
        self.end = start + 1
        self.type = 1
        self.seq = seq
        self.readnames = [read.query_name]
        self.uniqStr = "Insertion|%s:%s,%s|%s" % (self.contig, self.start,
                                                  self.end, self.seq)
        try:
            self.shen = min([shen(handle.fetch(read.reference_id,
                                               start=start - window,
                                               end=start)),
                             shen(handle.fetch(read.reference_id,
                                               start=self.end,
                                               end=self.end + window))])
        except AttributeError:
            self.shen = -1.
        self.handle = handle
        self.shenwindow = window
        self.StartStops = [ssStringFromRead(read)]

    @cython.returns(cython.str)
    def __str__(self):
        return self.uniqStr + "|%s" % len(self.readnames)



cdef class Deletion(AbstractIndelContainer):
    """
    Expansion of the AbstractIndelContainer object.
    "Start" is the first missing base.
    "End" is the last missing base.
    If the deletion is of length one, start and end should be the same.
    """

    def __init__(self, pysam.calignmentfile.AlignedSegment read,
                 cython.str contig=None, cython.long start=-1,
                 cython.long end=-1,
                 pysam.cfaidx.FastaFile handle=None, cython.long window=20):
        self.contig = contig
        self.start = start
        self.end = end
        self.type = -1
        self.readnames = [read.query_name]
        self.uniqStr = "Deletion|%s:%s,%s" % (self.contig, self.start,
                                              self.end)
        try:
            self.shen = min([shen(handle.fetch(contig,
                                               start=start - window,
                                               end=start)),
                             shen(handle.fetch(contig,
                                               start=end, end=end + window))])
        except AttributeError:
            self.shen = -1.
        self.handle = handle
        self.shenwindow = window
        self.StartStops = [ssStringFromRead(read)]

    @cython.returns(cython.str)
    def __str__(self):
        return self.uniqStr + "|%s" % len(self.readnames)


cdef class IndelQuiver(object):

    """
    Class for holding on to the indel objects. Holds a list for each
    type of indel, as well as a dictionary where the keys are unique
    string descriptors for the mutation and the values are lists
    of read names for reads supporting that mutation.
    Counts is a similar object, but with the length of the data
    field as a value instead of the list itself.
    """
    def __init__(self, cython.str ref=None, cython.long window=10,
                 cython.long minMQ=0, cython.long minFM=0,
                 cython.str bam=None, cython.long minNumSS=0,
                 cython.float minShen=minShen, cython.long minPairs=minPairs):
        self.data = {}
        self.readnames = {}
        self.counts = {}
        self.fastaRef = pysam.FastaFile(ref)
        self.bam = pysam.AlignmentFile(bam, "rb")
        self.window = window
        self.minMQ = minMQ
        self.minFM = minFM
        self.minPairs = minPairs
        self.minShen = minShen
        self.minNumSS = minNumSS


    @cython.returns(cython.long)
    def __len__(self):
        return len(self.data)

    @cython.returns(list)
    def __getitem__(self, cython.str key):
        return self.data[key]

    def __setitem__(self, cython.str key, list value):
        self.data[key] = value

    def setIndelShen(self, AbstractIndelContainer_t indelObj):
        indelObj.shen = min([shen(self.fastaRef(indelObj.contig,
                                                indelObj.start - self.window,
                                                indelObj.start)),
                             shen(self.fastaRef(indelObj.contig, indelObj.end,
                                                indelObj.end + self.window))])
        indelObj.shenwindow = self.window

    def iterkeys(self):
        return self.data.iterkeys()

    @cython.returns(list)
    def keys(self):
        return self.data.keys()

    def iteritems(self):
        return self.data.iteritems()

    @cython.returns(list)
    def items(self):
        return self.data.items()

    def itervalues(self):
        return self.data.itervalues()

    @cython.returns(list)
    def values(self):
        return self.data.values()

    def addRead(self, pysam.calignmentfile.AlignedSegment read):
        cdef cython.str SVTag
        cdef AbstractIndelContainer_t Indel
        if(read.opt("FM") < self.minFM):
            return
        if(read.mapping_quality < self.minMQ):
            return
        try:
            SVTag = read.opt("SV")
        except KeyError:
            raise ThisIsMadness("read must have an SV tag for IndelQuiver!")
        if("DSI" in SVTag):
            Indel = GetInsertionFromAlignedSegment(read, handle=self.fastaRef)
            try:
                self[Indel.uniqStr].register(read)
            except KeyError:
                self.setIndelShen(Indel)
                self[Indel.uniqStr] = Indel
        if("DSD" in SVTag):
            Indel = GetDeletionFromAlignedSegment(read, handle=self.fastaRef)
            try:
                self[Indel.uniqStr].register(read)
            except KeyError:
                self.setIndelShen(Indel)
                self[Indel.uniqStr] = Indel
        self.counts = {key: len(values) for key, values in
                       self.readnames.itervalues()}


    def mergeQuiver(self, IndelQuiver_t quiverObj):
        cdef cython.str key
        for key in quiverObj.readnames.keys():
            try:
                self.readnames[key] += quiverObj[key]
            except KeyError:
                self.readnames[key] = quiverObj[key]
        self.counts = {key: len(values) for key, values in
                       self.readnames.itervalues()}
        for key in quiverObj.keys():
            try:
                self[key].merge(quiverObj[key])
            except KeyError:
                self[key] = quiverObj[key]

    @cython.returns(dict)
    def getIndelCounter(self, cython.str uniqStr):
        try:
            return cyfreq(self.readnames[uniqStr])
        except KeyError:
            return {}

    def makeVCFLine(self, AbstractIndelContainer_t IC):
        return IDVCFLine(IC, self)


class IDVCFLine:

    """
    Contains data for writing to a VCF Line, where each ALT has its own line.

    minNumSS is the minumum number of start/stop combinations required to
    support a variant call.
    """
    def __init__(self, AbstractIndelContainer_t IC, IndelQuiver_t quiver=None):
        cdef pysam.calignmentfile.PileupColumn PileupCol
        cdef pysam.calignmentfile.IteratorColumnRegion pileupIt
        cdef cython.long tmpCov
        cdef cython.float MDP
        cdef tuple i
        self.CHROM = IC.contig
        self.NumStartStops = len(set(IC.StartStops))
        NDP = len(set(IC.readnames))
        if(isinstance(IC, Insertion)):
            """
            I would subtract one to get the start, but I would also
            add one to correct for 0 vs 1-based.
            One is subtracted from the fetch.
            """
            self.TYPE = "ins"
            self.POS = IC.start
            self.REF = quiver.fastaRef.fetch(IC.contig,
                                             self.POS - 1,
                                             self.POS)
            self.ALT = self.REF + IC.seq
            self.LEN = len(IC.seq)
        elif(isinstance(IC, Deletion)):
            """
            Because the start of a Deletion object is the first deleted
            base, the start of the variant is that minus one.
            Then, we're incrementing for 0 vs 1
            The end of the fetch is incremented by one for 0 vs 1 and
            once more to get the base after the last deleted base.
            """
            self.TYPE = "del"
            self.POS = IC.start
            self.REF = quiver.fastaRef.fetch(IC.contig,
                                             self.POS - 1,
                                             IC.end + 1)
            self.ALT = self.REF[0]
            self.LEN = IC.end - IC.start
        else:
            raise ThisIsMadness("Sorry, I haven't finished this "
                                "VCF writer for complex indels.")
        self.ID = IC.uniqStr
        self.reverseStrandFraction = sum(
            ["reverse" in ssString for ssString in
             IC.StartStops]) / len(IC.StartStops)
        self.BothStrandSupport = (self.reverseStrandFraction > 0.001 and
                                  self.reverseStrandFraction < 0.999)
        self.QUAL = 20 * len(IC)  # 20 for each supporting read. Will change!
        pileupIt = quiver.bam.pileup(self.CHROM, self.POS, max_depth=200000)
        PileupCol = pileupIt.next()
        while PileupCol.pos < self.POS - 5:
            PileupCol = pileupIt.next()
        tmpCov = 0
        while PileupCol.pos < self.POS + 5:
            tmpCov += PileupCol.n
            PileupCol = pileupIt.next()
        self.MDP = tmpCov * 1. / 10
        self.FILTER = None
        self.NDPS = len([i for i in cyfreq(IC.readnames).iteritems() if
                         i[1] > 1])
        self.DPA = len(IC.readnames)
        if(NDP < quiver.minPairs):
            try:
                self.FILTER.append(";InsufficientReadPairs")
            except AttributeError:
                self.FILTER = "InsufficientReadPairs"
        if(self.NumStartStops < quiver.minNumSS):
            try:
                self.FILTER.append(";InsufficientStopStarts")
            except:
                self.FILTER = "InsufficientStopStarts"
        if(IC.shen < quiver.minShen):
            try:
                self.FILTER.append(";LowComplexity")
            except:
                self.FILTER = "LowComplexity"
        if(self.FILTER is None):
            self.FILTER = "PASS"
        self.InfoFields = {"SHENWINDOW": quiver.window,
                           "MINSHEN": quiver.minShen}
        self.FormatFields = {"SHEN": IC.shen, "TYPE": self.TYPE,
                             "NDPS": self.NDPS, "DPA": self.DPA,
                             "MDP": self.MDP, "LEN": self.LEN}
        ffkeys = sorted(self.FormatFields.keys())
        self.FormatStr = (
            ":".join(ffkeys) +
            "\t" + ":".join(str(
                self.FormatFields[key]) for key in ffkeys))
        self.str = "\t".join(map(str, [self.CHROM,
                                       self.POS, self.ID,
                                       self.CONS, self.ALT,
                                       self.QUAL, self.FILTER,
                                       self.InfoStr, self.FormatStr]))

    def update(self):
        ffkeys = sorted(self.FormatFields.keys())
        self.FormatKey = ":".join(ffkeys)
        self.FormatValue = ":".join([str(self.FormatFields[key])
                                     for key in
                                     ffkeys])
        self.FormatStr = self.FormatKey + "\t" + self.FormatValue
        self.InfoStr = ";".join([key + "=" + str(self.InfoFields[key])
                                 for key in sorted(self.InfoFields.keys())])

    def __str__(self):
        self.update()
        self.str = "\t".join(map(str, [self.CHROM, self.POS,
                                       self.ID, self.REF, self.ALT,
                                       self.QUAL, self.FILTER, self.InfoStr,
                                       self.FormatStr]))
        return self.str

# Info fields
IDInfoDict = {}
IDInfoDict["MINSHEN"] = HeaderInfoLine(
    ID="MINSHEN", Type="Float", Number="1",
    Description=("Minimum Shannon entropy for both"
                 " preceding and succeeding regions"))
IDInfoDict["SHENWINDOW"] = HeaderInfoLine(
    ID="SHENWINDOW", Type="Integer", Number="1",
    Description=("Number of preceding or succeeding bases to include in "
                 "Shannon entropy calculations."))


# Filter fields
IDFilterDict = {}
IDFilterDict["PASS"] = HeaderFilterLine(ID="PASS",
                                        Description="All filters passed")
IDFilterDict["LowComplexity"] = HeaderFilterLine(
    ID="LowComplexity",
    Description=("Variant's flanking regions have Shannon entropy "
                 "below a required threshold."))
IDFilterDict["InsufficientReadPairs"] = HeaderFilterLine(
    ID="InsufficientReadPairs",
    Description=("Variant not supported by sufficient "
                 "concordant pairs of duplex reads"))
IDFilterDict["InsufficientStopStarts"] = HeaderFilterLine(
    ID="InsufficientStopStarts",
    Description="Number of start/stop coordinates is less than minNumSS")

# Format fields
IDFormatDict = {}
IDFormatDict["TYPE"] = HeaderFormatLine(
    Type="String", ID="TYPE",
    Description=("The type of allele, either snp, mnp, "
                 "ins, del, tra, or complex."),
    Number="A")
IDFormatDict["NDPS"] = HeaderFormatLine(
    Type="Integer", ID="NDPS", Number="A",
    Description=("Number of duplex reads (where read 1 and read 2 map to the "
                 "same location) supporting variant from both directions."))
IDFormatDict["DPA"] = HeaderFormatLine(
    Type="Integer", ID="AC", Number="A",
    Description=("Number of reads passing filters supporting variant."))
IDFormatDict["SHEN"] = HeaderFormatLine(
    ID="SHEN", Type="Float", Number="A",
    Description=("Shannon entropy of preceding or succeeding "
                 "sequence, whichever is smaller."))
IDFormatDict["LEN"] = HeaderFormatLine(
    ID="LEN", Description="Length of allele (deletion or insertion)",
    Type="Integer", Number="A")
IDFormatDict["MDP"] = HeaderFormatLine(
    ID="MDP", Number="A", Type="Float",
    Description=("Mean depth of coverage for region in which indel was "
                 "called. DOC is calculated by averaging 5 bases before and"
                 " 5 bases after the event to give context for the event."))


def GetIDVCFHeader(fileFormat="default", commandStr="default",
                   reference="default",
                   header="default", sampleName="DefaultSampleName"):
    reference = reference.split("/")[-1]
    # If the reference is a path, trim to just the name of the file.
    HeaderLinesStr = ""
    # fileformat line
    HeaderLinesStr += str(HeaderFileFormatLine(
        fileformat=fileFormat)) + "\n"
    # FILTER lines
    for key in sorted(IDFilterDict.keys()):
        HeaderLinesStr += str(IDFilterDict[key]) + "\n"
    # INFO lines
    HeaderLinesStr += "\n".join([str(IDInfoDict[key]) for
                       key in sorted(IDInfoDict.keys())]) + "\n"
    # FORMAT lines
    HeaderLinesStr += "\n".join([str(IDFormatDict[key]) for
                       key in IDFormatDict.keys()]) + "n"
    # commandline line
    if(commandStr != "default"):
        HeaderLinesStr += str(HeaderCommandLine(
            commandStr=commandStr)) + "\n"
    # reference line
    if(reference != "default"):
        HeaderLinesStr += str(HeaderReferenceLine(
            reference=reference)) + "\n"
    # contig lines
    if(header != "default"):
        HeaderLinesStr += GetContigHeaderLines(header) + "\n"
    HeaderLinesStr += "\t".join(["#CHROM", "POS",
                                 "ID", "REF",
                                 "ALT", "QUAL",
                                 "FILTER", "INFO", "FORMAT",
                                 sampleName]) + "\n"
    return HeaderLinesStr