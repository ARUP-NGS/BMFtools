# cython: c_string_type=str, c_string_encoding=ascii
## Standard Library module-level imports
from __future__ import division
import logging
import sys
import uuid

## Third party imports
import pysam
import cython
import numpy as np

## Standard Library local imports
from array import array
from itertools import izip, groupby
from operator import attrgetter as oag, methodcaller as mc
from subprocess import check_call
from sys import maxint

##BMFTools imports
from utilBMF.ErrorHandling import ThisIsMadness
from utilBMF.HTSUtils import printlog as pl, BamTag, TrimExt

##DEFINES
oagsk = oag("firstMapped")
omcfp = mc("cGetRefPosForFirstPos")
oagbase = oag("base")
oagop = oag("operation")
oagqual = oag("quality")
oagag = oag("agreement")
oagtag = oag("tag")
chrDict = {x: chr(x) for x in xrange(126)}


"""
Contains utilities for merging a pair of overlapping alignments
into a single read.
"""

## DEFINES
# Maps CigarOp numbers to ASCII characters
# and the reverse
CigarDict = {0: 77, 1: 73, 2: 68, 3: 83, 4: 78, 5: 72,
             6: 80, 7: 61, 8: 88,
             77: 0, 73: 1, 68:2, 83: 3, 78: 4, 72: 5,
             80: 6, 61: 7, 88: 8}

cdef class ListBool:
    """
    Used to strongly type return type of a list and a bool
    """
    def __cinit__(self, list List, bint Bool):
        self.List = List
        self.Bool = Bool
        


cdef list CigarOpToLayoutPosList(int offset, int cigarOp, int cigarLen,
                                 pysam.calignmentfile.AlignedSegment rec):
    cdef object x0, x1
    cdef char CigarChar
    cdef ndarray[long, ndim=1] quals, agrees
    '''
    First case - 'M'
    Second case - 'I'
    Third case - 'S'
    Fourth case - 'D'
    '''
    CigarChar = CigarDict[cigarOp]
    try:
        quals = np.array(rec.opt("PV").split(","), dtype=np.int64)
    except KeyError:
        pl("Watch out - PV tag not set.", level=logging.DEBUG)
        quals = np.array(rec.query_qualities, dtype=np.int64)
        # Let's make sure that these don't need reversal, too!
    try:
        agrees = np.array(rec.opt("FA").split(","), dtype=np.int64)
    except KeyError:
        pl("Watch out - FA tag not set.", level=logging.DEBUG)
        agrees = np.array([1] * len(rec.sequence), dtype=np.int64)
    return [LayoutPos(pos=x1, readPos=x0, operation=CigarChar,
                      base=ord(rec.seq[x0]), quality=quals[x0],
                      agreement=agrees[x0]) if
            x1 is not None and x0 is not None else
            LayoutPos(-1, x0, operation=CigarChar,
                      base=ord(rec.seq[x0]),
                      quality=quals[x0], agreement=agrees[x0]) if
            x0 is not None else
            LayoutPos(x1, x0, operation=83,  # 83 == "S"
                      base=ord(rec.seq[x0]),
                      quality=quals[x0], agreement=agrees[x0]) if
            CigarChar == 83 else
            LayoutPos(x1, -1, operation=68,
                      base=68,  # Base set to "D"
                      quality=-1, agreement=-1)
            for x0, x1 in rec.aligned_pairs[offset:offset + cigarLen]]


@cython.returns(Layout_t)
def makeLayout(pysam.calignmentfile.AlignedSegment rec):
    return Layout(makeLayoutTuple)


@cython.returns(tuple)
def makeLayoutTuple(pysam.calignmentfile.AlignedSegment rec):
    cdef tuple pair
    cdef int offset, firstMapped, i
    cdef list PosList
    offset = 0
    firstMapped = -1
    PosList = []
    for i, pair in enumerate(rec.cigar):
        PosList += CigarOpToLayoutPosList(offset, pair[0], pair[1], rec)
        if(pair[0] == 0 and firstMapped < 0):  # 0->"M"
            firstMapped = offset  # Facilitates coordinating merging pairs.
        offset += pair[1]
    return (rec, PosList, firstMapped)


@cython.returns(bint)
def lambda1None(tuple i):
    """
    If a tuple in a cigar returns true for this function,
    then that base is either deleted or soft-clipped, which
    you can tell based on whether it is in the middle or the end of a read.
    """
    return i[1] is None


@cython.returns(int)
def getFirstMappedRefPos(pysam.calignmentfile.AlignedSegment rec):
    cdef tuple i
    return [i for i in rec.aligned_pairs if i[1] is not None][0][1]


cdef class LayoutPos:
    """
    Holds one layout position - either part of the reference,
    part of the read, or both.

    All of these fields are chars/ints.
    pos:
        -1 if "I" or "S"
    readPos:
        -1 if "D"
    mergeAgreed:
        -1 = Not merged
        0 = Merge did not agree
        1 = Merge Agreed
    base:
        ascii char/int for the base
        [ord(i) for i in 'ACGNT']
        [65, 67, 71, 78, 84]
    operation:
        ascii char/int for the operation
        [ord(i) for i in "MIDNSHP=X"]
        [77, 73, 68, 78, 83, 72, 80, 61, 88]
    """
    def __init__(self, int pos=-1, int readPos=-1,
                 char base=-1, char operation=-1,
                 int quality=-1, int agreement=-1,
                 bint isMerged=False, char mergeAgreed=-1):
        self.pos = pos
        self.readPos = readPos
        self.operation = operation
        self.base = base
        self.quality = quality if(self.base != 78) else 0  # 78 == "N"
        self.agreement = agreement
        self.isMerged = isMerged
        self.mergeAgreed = mergeAgreed

    cpdef bint ismapped(self):
        return self.operation == 77  # 77 == "M"

    cdef bint getMergeAgreed(self):
        return self.mergeAgreed > 0

    @cython.returns(cystr)
    def __str__(self):
        return "%s|%s|%s|%s|%s|%s|%s|%s" % (
            self.pos, self.readPos, chr(self.base), chr(self.operation),
            self.quality, self.agreement, self.isMerged, self.mergeAgreed)

cdef class ArrayLayoutPos:
    # pos, readPos, quality, agreement, operation, base
    # cdef public cystr operation, base

    def __cinit__(self, int pos=-1, int readPos=-1,
                  int base=-1, int operation=-1,
                  int quality=-1, int agreement=-1):
        self.values = array('i', [pos, readPos, base,
                                  operation, quality, agreement])
        raise NotImplementedError

    def __init__(self, *args):
        self.__cinit__(*args)

    def __getitem__(self, int index):
        return self.values[index]

    cpdef bint ismapped(self):
        return self[3] == 77

    @cython.returns(cystr)
    def __str__(self):
        return "%s|%s|%s|%s|%s|%s" % (self.values[0], self.values[1],
                                      chr(self.values[2]),
                                      chr(self.values[3]),
                                      self.values[4], self.values[5])


cdef class Layout(object):
    """
    Holds a read and its layout information.

    This doctest was written so that it would load in one read,
    make the string, and then hash that value. Since it wouldn't
    match the interpreter's output to have a gigantic line and it would have
    to violate pep8, I decided to test the value by its hash rather than by
    string agreement.
    Unfortunately, this only works on 64-bit systems, as the hash function
    is different for 32-bit systems.
    >>> from sys import maxint
    >>> from pysam import AlignmentFile as af
    >>> handle = af("utilBMF/example.bam", "rb")
    >>> returnStr = str(Layout.fromread(handle.next()))
    >>> hashreturn = -8225309399721982299
    >>> hash(returnStr)
    -8225309399721982299
    """

    @classmethod
    def fromread(cls, pysam.calignmentfile.AlignedSegment rec):
        return cls(*makeLayoutTuple(rec))

    def __getitem__(self, index):
        return self.positions[index]

    def __len__(self):
        return len(self.positions)

    cdef ndarray[char] getSeqArr(self, dict chrDict=chrDict):
        """Returns a character array of the base calls
        if the base calls aren't "S" (83) or "D" (68)
        """
        cdef char i
        cdef LayoutPos_t pos
        return np.char.array([chrDict[pos.base]
                              for pos in self.positions if
                              pos.operation != 83 and pos.operation != 68 and
                              pos.agreement > -1],
                             itemsize=1)

    cpdef cystr getSeq(self):
        return self.getSeqArr().tostring()

    cdef int cGetRefPosForFirstPos(self):
        """cdef class wrapped by pGetRefPosForFirstPos
        """
        cdef LayoutPos_t i
        cdef int count
        for count, i in enumerate(self):
            if(i.operation == "M"):
                return i.pos - count

    cpdef int pGetRefPosForFirstPos(self):
        return self.cGetRefPosForFirstPos()

    cpdef int getAlignmentStart(self):
        cdef LayoutPos_t i
        for i in self.positions:
            if(i.operation == 77):  # operation is "M"
                return i.pos

    cpdef ndarray[int, ndim=1] getAgreement(self):
        """cpdef wrapper of cGetAgreement
        """
        return self.cGetAgreement()

    cdef ndarray[int, ndim=1] cGetAgreement(self):
        cdef int i
        cdef LayoutPos_t pos
        # Ask if >= 0. My tests say it's ~1% faster to ask (> -1) than (>= 0).
        return np.array([i for
                         i in [pos.agreement for
                               pos in self.positions] if i > -1],
                        dtype=np.int64)

    cdef ndarray[int, ndim=1] cGetQual(self):
        cdef LayoutPos_t pos
        # Ask if >= 0. My tests say it's ~1% faster to ask (> -1) than (>= 0).
        # pos.base == 66 for "B", which is a blank spot.
        # quality is set to less than 0 for an "N" cigar operation.
        return np.array([pos.quality for
                         pos in self.positions if pos.base != 66 and
                         pos.quality > -1],
                        dtype=np.int64)

    cpdef ndarray[int, ndim=1] getQual(self):
        return self.cGetQual()

    cdef cystr cGetQualString(self, dict ph2chrDict=ph2chrDict):
        cdef int i
        return "".join([ph2chrDict[i] for i in self.getQual()])

    cpdef getQualString(self):
        return self.cGetQualString()

    cdef int cGetLastRefPos(self):
        """cdef version of getLastRefPos
        """
        cdef LayoutPos_t i
        cdef int lastMPos, countFromMPos, count
        lastMPos = 0
        for count, i in enumerate(self):
            countFromMPos += 1
            if(i.operation == 77):  # 77 == M in ASCII
                lastMPos = i.pos
                countFromMPos = 0
        return lastMPos + countFromMPos

    cpdef int getLastRefPos(self):
        """
        Finds the reference position which "matches" the last base in the
        layout. This counts for any cigar operation.
        """
        return self.cGetLastRefPos()

    cpdef ndarray[char, ndim=1] getMergedPositions(self):
        """
        Returns the indices of the bases in the merged read
        which have been merged.
        """
        return self.cGetMergedPositions()

    cpdef ndarray[char, ndim=1] getMergeAgreements(self):
        """
        Returns the indices of the bases in the merged read
        which have been merged successfully.
        """
        return self.cGetMergeAgreements()

    cdef ndarray[char, ndim=1] cGetMergedPositions(self):
        """
        Returns the indices of the bases in the merged read
        which have been merged.
        """
        cdef LayoutPos_t lpos
        return np.array([lpos.readPos for lpos in self.positions if
                         lpos.isMerged],
                        dtype=np.int8)

    cdef ndarray[char, ndim=1] cGetMergeAgreements(self):
        """
        Returns the indices of the bases in the merged read
        which have been merged successfully.
        """
        cdef LayoutPos_t lpos
        return np.array([lpos.readPos for lpos in self.positions if
                         lpos.getMergeAgreed()],
                        dtype=np.int8)

    cdef ndarray[int, ndim=1] cGetGenomicDiscordantPositions(self):
        """
        Returns a 1-d array of integers for genomic positions which
        were discordant between the pairs.
        """
        cdef LayoutPos_t p
        if(self.isMerged is False):
            return np.array([], dtype=np.int64)
        return np.array([p.pos for p in self.positions if
                         p.mergeAgreed == 0 and p.pos >= 0],
                        dtype=np.int64)

    cdef ndarray[int, ndim=1] cGetReadDiscordantPositions(self):
        """
        Returns a 1-d array of integers for genomic positions which
        were discordant between the pairs.
        """
        cdef LayoutPos_t p
        if(self.isMerged is False):
            return np.array([], dtype=np.int64)
        return np.array([p.readPos for p in self.positions if
                         p.mergeAgreed == 0 and p.readPos >= 0],
                        dtype=np.int64)

    cdef update_tags_(self):
        self.tagDict["PV"] = BamTag(
            "PV", "Z", ",".join(self.getQual().astype(str)))
        self.tagDict["FA"] = BamTag(
            "FA", "Z", ",".join(self.getAgreement().astype(str)))
        if(self.isMerged):
            self.tagDict["PM"] = BamTag(
                "PM", "Z", ",".join(self.getMergedPositions().astype(str)))
            self.tagDict["MA"] = BamTag(
                "MA", "Z", ",".join(self.getMergeAgreements().astype(str)))
            self.tagDict["DG"] = BamTag(
                "DG", "Z", ",".join(
                    self.cGetGenomicDiscordantPositions().astype(str)))
            self.tagDict["DR"] = BamTag(
                "DR", "Z", ",".join(
                    self.cGetReadDiscordantPositions().astype(str)))
            # Update it for the merged world!
            # Original template length
            self.tagDict["ot"] = BamTag("ot", "i", self.tlen)
            # Original mate position
            self.tagDict["mp"] = BamTag("mp", "i", self.pnext)
            # Original mapping quality
            self.tagDict["om"] = BamTag("om", "i", self.mapq)
            # Original mapped position
            self.tagDict["op"] = BamTag("op", "i", self.InitPos)
            self.tagDict["MP"] = BamTag("MP", "A", "T")

    cpdef update_tags(self):
        self.update_tags_()

    def update(self):
        cdef LayoutPos_t pos
        cdef int count
        self.update_tags()
        if(self.isMerged):
            self.tlen = len(self.getSeqArr())
            self.pnext = 0
            self.mapq = -1
            self.rnext = "*"
            self.flag = 2 + (16 if(self.is_reverse) else 32)
            for count, pos in enumerate(self):
                pos.readPos = count

    @cython.returns(ndarray)
    def getOperations(self, oagop=oagop):
        """
        In [15]: %timeit [chr(p.operation) for p in l1.positions]
        The slowest run took 5.12 times longer than the fastest. This could mean that an intermediate result is being cached 
        100000 loops, best of 3: 9.35 us per loop

        In [16]: %timeit map(chr, [p.operation for p in l1.positions])
        The slowest run took 6.09 times longer than the fastest. This could mean that an intermediate result is being cached 
        100000 loops, best of 3: 8.69 us per loop
        """
        cdef LayoutPos_t pos
        # Skip this operation if the operation is "N"
        return np.array(map(chr, [pos.operation for pos in
                                  self.positions if pos.operation != 78]))

    cdef cystr cGetCigarString(self):
        return "".join([str(len(list(g))) + k for
                        k, g in groupby(self.getOperations())])

    cpdef cystr getCigarString(self):
        return self.cGetCigarString()

    @cython.returns(list)
    def get_tags(self, oagtag=oagtag):
        self.update_tags()
        return sorted(self.tagDict.itervalues(), key=oagtag)

    def getFlag(self):
        self.update()
        return self.flag

    @cython.returns(cystr)
    def __str__(self):
        """
        Converts the record into a SAM record.
        Note: the position is incremented by 1 because SAM positions are
        1-based instead of 0-based.
        """
        self.update()
        return "\t".join(map(
                str, [self.Name, self.getFlag(), self.contig,
                      self.getAlignmentStart() + 1, self.mapq,
                      self.getCigarString(), self.rnext, self.pnext + 1,
                      self.tlen, self.getSeq(), self.getQualString()] +
                self.get_tags()))

    def __init__(self, pysam.calignmentfile.AlignedSegment rec,
                 list layoutPositions, int firstMapped,
                 PysamToChrDict=PysamToChrDict):
        cdef tuple tag
        self.mapq = rec.mapq
        self.read = rec
        self.positions = layoutPositions
        self.firstMapped = firstMapped
        self.InitPos = rec.pos
        self.Name = rec.query_name
        self.contig = PysamToChrDict[rec.reference_id]
        self.flag = rec.flag
        # When the get_tags(with_value_type=True) is implemented,
        # then switch the code over.
        # Then I can change the way I make BAM tags, since
        # with_value_type is the optional argument to add to get_tags
        # on a pysam.calignmentfile.AlignedSegment object.
        # This will ensure that I don't need to already know
        # the type beforehand. (IE, get rid of the type dict)
        self.tagDict = {tag[0]: BamTag.fromtuple(tag) for tag
                        in rec.get_tags() if tag[0] not in ["PV", "FA"]}
        self.rnext = PysamToChrDict[rec.mrnm]
        self.pnext = rec.mpos
        self.tlen = rec.tlen
        self.isMerged = (rec.has_tag("MP") and rec.opt("MP") == "T")
        self.is_reverse = rec.is_reverse


cdef LayoutPos_t cMergePositions(LayoutPos_t pos1, LayoutPos_t pos2):
    """Merges two positions. Order does matter - pos1 overrides pos2 when
    pos2 is soft-clipped.
    """
    if(pos1.operation != pos2.operation):
        if(pos2.operation == 83):  # if pos2.operation is "S"
            if(pos1.base == pos2.base):
                return LayoutPos(pos1.pos, pos1.readPos, pos1.base,
                                 pos1.operation,
                                 pos1.quality + pos2.quality,
                                 pos1.agreement + pos2.agreement,
                                 isMerged=True, mergeAgreed=1)
            else:
                return LayoutPos(pos1.pos, pos1.readPos, 66, 78, -137, -137,
                                 isMerged=True, mergeAgreed=0)
                # B marks a "Blank" space - for these bases, no quality entry
                # and no sequence entry are added to the output sam read.
                # N is the Cigar Operation
                # Quality 0, agreement 0.
    if(pos1.base == pos2.base):
        if(pos1.operation == pos2.operation):
            return LayoutPos(pos1.pos, pos1.readPos, pos1.base, pos1.operation,
                             pos1.quality + pos2.quality,
                             pos1.agreement + pos2.agreement, isMerged=True,
                             mergeAgreed=1)
        elif(pos2.operation == 83):  # if pos2.operation is "S"
            return LayoutPos(pos1.pos, pos1.readPos, pos1.base,
                             pos1.operation,
                             pos1.quality + pos2.quality,
                             pos1.agreement + pos2.agreement,
                             isMerged=True, mergeAgreed=1)
        else:
            return LayoutPos(pos1.pos, pos1.readPos, 66, 78, -137, -137,
                             isMerged=True, mergeAgreed=0)
    else:
        if(pos1.operation == pos2.operation):
            if(pos1.quality > pos2.quality):
                return LayoutPos(
                    pos1.pos, pos1.readPos, pos1.base, pos1.operation,
                    pos1.quality - pos2.quality, pos1.agreement,
                    isMerged=True, mergeAgreed=0)
            else:
                return LayoutPos(
                    pos1.pos, pos1.readPos, pos2.base, pos1.operation,
                    pos2.quality - pos1.quality, pos2.agreement,
                    isMerged=True, mergeAgreed=0)
        else:
            return LayoutPos(pos1.pos, pos1.readPos, 66, 78, -137, -137,
                             isMerged=True, mergeAgreed=0)


cdef ListBool cMergeLayoutsToList(Layout_t L1, Layout_t L2):
    """
    Merges two Layouts into a list of layout positions.

    First, it takes the positions before the overlap starts.
    Then it merges the overlap position by position.
    Then it takes the positions after the overlap ends.
    :param Layout_t L1: One Layout object
    :param Layout_t L2: Another Layout object
    :param oagsk: A method caller function. Locally redefined for speed.

    :return list Merged Positions
    :return bool Whether the merge was successful
    """
    cdef int offset
    cdef Layout_t tmpPos
    cdef LayoutPos_t pos1, pos2
    if(LayoutsOverlap(L1, L2) is False):
        return ListBool([], False)
    if(L1.cGetRefPosForFirstPos() > L2.cGetRefPosForFirstPos()):
        tmpPos = L1
        L1 = L2
        del tmpPos
    # L1, L2 = sorted((L1, L2), key=omcfp) previous python code
    # Rewritten to avoid the python object omcfp
    offset = L2.cGetRefPosForFirstPos() - L1.cGetRefPosForFirstPos()
    try:
        return ListBool((L1[:offset] + [cMergePositions(pos1, pos2) for
                                        pos1, pos2 in izip(L1[offset:], L2)] +
                         L2[len(L1) - offset:]), True)
        '''
        return (L1[:offset] +
                [MergePositions(pos1, pos2) for
                 pos1, pos2 in izip(L1[offset:], L2)] +
                L2[len(L1) - offset:]), True
        '''
    except ThisIsMadness:
        pl("ThisIsMadness got thrown in cMergeLayoutsToList!"
            "Layout 1: %s. Layout 2: %s" % (str(L1), str(L2)),
            level=logging.DEBUG)
        return ListBool(L1[:] + L2[len(L1) - offset:], False)


cpdef Layout_t MergeLayoutsToLayout(Layout_t L1, Layout_t L2):
    cdef list layoutList
    cdef cystr Name
    cdef bint Success
    cdef ListBool ret
    ret = cMergeLayoutsToList(L1, L2)
    if(ret.Bool is False):
        L1.tagDict["MP"] = BamTag("MP", tagtype="Z", value="F")
        L2.tagDict["MP"] = BamTag("MP", tagtype="Z", value="F")
        return None
    L1.positions = ret.List
    L1.tagDict["MP"] = BamTag("MP", tagtype="Z", value="T")
    L1.isMerged = True
    L1.update()
    return L1


cpdef bint LayoutsOverlap(Layout_t L1, Layout_t L2):
    return cReadsOverlap(L1.read, L2.read)


def MergePairedAlignmentsToBam(cystr inBAM, cystr outBAM=None,
                               bint coorsort=False):
    cdef cystr uuidStr = str(uuid.uuid4(
        ).get_hex().upper()[0:8])
    if(outBAM is None):
        outBAM = TrimExt(inBAM) + ".PairMergeProcessed.bam"
    if(coorsort):
        cStr = ("python -c 'from utilBMF.MergePairedAlignments import MergePa"
                "iredAlignments;MergePairedAlignments(\"%s\", Pipe" % inBAM +
                "ToBam=True)' | samtools sort -T"
                " %s -o %s -O bam -" % (uuidStr, outBAM))
    else:
        cStr = ("python -c 'from utilBMF.MergePairedAlignments import Merge"
                "PairedAlignments;MergePairedAlignments(\"%s\", " % inBAM +
                "PipeToBam=True)' | samtools view -Sbh - > %s" % outBAM)
    pl("About to call MergePairedAlignmentsToBam. cStr: %s" % cStr)
    check_call(cStr, shell=True)
    return outBAM


def MergePairedAlignments(cystr inBAM, cystr outSAM=None,
                          bint PipeToBam=False):
    cdef cAlignedSegment read, read1, read2
    cdef Layout_t Layout1, Layout2, retLayout
    cdef int count = 0
    if(outSAM is None):
        outSAM = TrimExt(inBAM) + ".PairMergeProcessed.sam"
    inHandle = pysam.calignmentfile.AlignmentFile(inBAM, "rb")
    if(PipeToBam):
        outHandle = sys.stdout
    else:
        outHandle = pysam.calignmentfile.AlignmentFile(outSAM, "w",
                                                   template=inHandle)
    outStringHandle = open(outSAM, "w+")
    oshw = outStringHandle.write
    ohw = outHandle.write
    for read in inHandle:
        count += 1
        if(read.is_supplementary or read.is_secondary):
            ohw(read)
            continue
        if(read.is_read1):
            read1 = read
            continue
        read2 = read
        try:
            assert(read1.query_name == read2.query_name)
        except AssertionError:
            raise ThisIsMadness("Bam is either not name sorted or you are "
                                "missing a read from the pair around read "
                                "# %s in the bam." % count)
        Layout1 = Layout.fromread(read1)
        Layout2 = Layout.fromread(read2)
        retLayout = MergeLayoutsToLayout(Layout1, Layout2)
        if(retLayout is None):
            read1.setTag("MP", "Z", "F")
            read2.setTag("MP", "Z", "F")
            ohw(read1)
            ohw(read2)
            continue
        oshw(str(retLayout) + "\n")
    outStringHandle.close()
    outHandle.close()
    inHandle.close()
    return outSAM
