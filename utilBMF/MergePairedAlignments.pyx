# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True
from __future__ import division
import pysam
import cython
import numpy as np
import logging
from utilBMF.HTSUtils import printlog as pl, BamTag
from utilBMF.ErrorHandling import ThisIsMadness
from operator import attrgetter as oag, methodcaller as mc
from itertools import izip, groupby
from sys import maxint
from array import array
from utilBMF.HTSUtils cimport cReadsOverlap
oagsk = oag("firstMapped")
omcfp = mc("getRefPosForFirstPos")
oagbase = oag("base")
oagop = oag("operation")
oagqual = oag("quality")
oagag = oag("agreement")
oagtag = oag("tag")

chrDict = {x: chr(x) for x in xrange(126)}

cimport pysam.calignmentfile
cimport cython
cimport numpy as np


"""
Contains utilities for merging a pair of overlapping alignments
into a single read.
"""


@cython.returns(list)
def CigarOpToLayoutPosList(int offset, tuple cigarOp,
                           pysam.calignmentfile.AlignedSegment rec):
    cdef tuple x
    cdef cystr CigarChar
    cdef ndarray[long, ndim=1] quals, agrees
    '''
    First case - 'M'
    Second case - 'I'
    Third case - 'D'
    Fourth case - 'S'
    '''
    CigarChar = CigarDict[cigarOp[0]]
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
    return [LayoutPos(pos=x[1], readPos=x[0], operation=ord(CigarChar),
                      base=ord(rec.seq[x[0]]), quality=quals[x[0]],
                      agreement=agrees[x[0]]) if
            x[1] is not None and x[0] is not None else
            LayoutPos(-1, x[0], operation=ord(CigarChar),
                      base=ord(rec.seq[x[0]]),
                      quality=quals[x[0]], agreement=agrees[x[0]]) if
            x[0] is not None else
            LayoutPos(x[1], -1, operation=ord(CigarChar),
                      base=ord(CigarChar),
                      quality=quals[x[0]], agreement=agrees[x[0]]) if
            CigarChar == "S" else
            LayoutPos(x[1], -1, operation=ord(CigarChar),
                      base=ord(CigarChar),
                      quality=-1, agreement=-1)
            for x in rec.aligned_pairs[offset:offset + cigarOp[1]]]


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
        PosList += CigarOpToLayoutPosList(offset, pair, rec)
        if(CigarDict[pair[0]] == "M" and firstMapped < 0):
            firstMapped = offset  # Facilitates coordinating merging pairs.
        offset += pair[1]
    return (rec, PosList, firstMapped)


@cython.returns(cython.bint)
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

    """
    def __init__(self, int pos=-1, int readPos=-1,
                 int base=-1, int operation=-1,
                 int quality=-1, int agreement=-1):
        self.pos = pos
        self.readPos = readPos
        self.operation = operation
        self.base = base
        self.quality = quality if(self.base != 78) else 0  # 78 == "N"
        self.agreement = agreement

    cpdef cython.bint ismapped(self):
        return self.operation == 77  # 77 == "M"

    @cython.returns(cystr)
    def __str__(self):
        return "%s|%s|%s|%s|%s|%s" % (self.pos, self.readPos, chr(self.base),
                                      chr(self.operation), self.quality,
                                      self.agreement)


cdef class ArrayLayoutPos:
    # pos, readPos, quality, agreement, operation, base
    # cdef public cystr operation, base

    def __cinit__(self, int pos=-1, int readPos=-1,
                  int base=-1, int operation=-1,
                  int quality=-1, int agreement=-1):
        self.values = array('i', [pos, readPos, base,
                                  operation, quality, agreement])

    def __init__(self, *args):
        self.__cinit__(*args)

    def __getitem__(self, int index):
        return self.values[index]

    cpdef cython.bint ismapped(self):
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
        return np.char.array([chrDict[i] for i in
                              [pos.base for pos in
                               self.positions] if
                              i != 83 and i != 68])
    '''
    cdef ndarray[char] getSeqArr(self):
        """Returns a character array of the base calls
        if the base calls aren't "S" (83) or "D" (68)
        """
        cdef char i
        cdef LayoutPos_t pos
        return np.char.array(map(chr, [i for i in
                                       [pos.base for pos in
                                        self.positions] if
                                       i != 83 and i != 68]))
    '''

    cpdef cystr getSeq(self):
        return self.getSeqArr().tostring()

    cdef int getRefPosForFirstPos_(self):
        """cdef class wrapped by getRefPosForFirstPos
        """
        cdef LayoutPos_t i
        cdef int count
        for count, i in enumerate(self):
            if(i.operation == "M"):
                return i.pos - count

    cpdef int getRefPosForFirstPos(self):
        return self.getRefPosForFirstPos_()

    cpdef int getAlignmentStart(self):
        cdef LayoutPos_t i
        for i in self:
            if(i.operation == "M"):
                return i.pos

    cpdef ndarray[int, ndim=1] getAgreement(self):
        """cpdef wrapper of getAgreement_
        """
        return self.getAgreement_()

    cdef ndarray[int, ndim=1] getAgreement_(self):
        cdef int i
        cdef LayoutPos_t pos
        # Ask if >= 0. My tests say it's ~1% faster to ask (> -1) than (>= 0).
        return np.array([i for
                         i in [pos.agreement for
                               pos in self.positions] if i > -1],
                        dtype=np.int64)

    cdef ndarray[int, ndim=1] getQual_(self):
        cdef int i
        cdef LayoutPos_t pos
        # Ask if >= 0. My tests say it's ~1% faster to ask (> -1) than (>= 0).
        return np.array([i for
                         i in [pos.quality for
                               pos in self.positions] if i > -1],
                        dtype=np.int64)

    cpdef ndarray[int, ndim=1] getQual(self):
        return self.getQual_()

    cdef cystr getQualString_(self, dict ph2chrDict=ph2chrDict):
        cdef int i
        return "".join([ph2chrDict[i] for i in self.getQual()])

    cpdef getQualString(self):
        return self.getQualString_()

    cdef int getLastRefPos_(self):
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
        return self.getLastRefPos_()

    cdef update_tags_(self):
        self.tagDict["PV"] = BamTag(
            "PV", "Z", ",".join(self.getQual().astype(str)))
        self.tagDict["FA"] = BamTag(
            "FA", "Z", ",".join(self.getAgreement().astype(str)))
        self.tagDict["PM"] = BamTag(
            "PM", "Z", ",".join(self.getMergedPositions().astype(str)))

    cpdef update_tags(self):
        self.update_tags_()

    def update(self):
        cdef LayoutPos_t pos
        cdef int count
        if(self.isMerged):
            # Update it for the merged world!
            # Original template length
            self.tagDict["ot"] = BamTag("ot", "i", self.tlen)
            # Original mate position
            self.tagDict["mp"] = BamTag("mp", "i", self.pnext)
            # Original mapping quality
            self.tagDict["om"] = BamTag("om", "i", self.mapq)
            # Original mapped position
            self.tagDict["op"] = BamTag("op", "i", self.InitPos)
            self.tagDict["FM"].value *= 2  # Double the FM tag.
            self.tlen = 0
            self.pnext = 0
            self.mapq = -1
            self.tagDict["MP"] = BamTag("MP", "A", "T")
            self.rnext = "*"
            self.flag = 2 + (16 if(self.is_reverse) else 32)
            for count, pos in enumerate(self):
                pos.readPos = count
        self.update_tags()

    @cython.returns(ndarray)
    def getOperations(self, oagop=oagop):
        cdef LayoutPos_t pos
        return np.array(map(chr, [pos.operation for pos in self.positions]))

    cdef cystr getCigarString_(self):
        return "".join([str(len(list(g))) + k for
                        k, g in groupby(self.getOperations())])

    cpdef cystr getCigarString(self):
        return self.getCigarString_()

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


def LayoutSortKeySK(x, oagsk=oagsk):
    """
    Pre-allocating the sort key function is 20% as fast as allocating
    it at each call. Total sort is ~10% faster by preallocating.

    Faster yet, though, is calling oagsk directly, which is a little over
    an additional 20% speed increase for its use with map.
    """
    return oagsk(x)
'''
cdef LayoutPos_t MergePositions(LayoutPos p1, p2):
'''

cdef LayoutPos_t MergePositions(LayoutPos pos1, LayoutPos pos2):
    """Merges two positions. Order does matter - pos1 overrides pos2 when
    pos2 is soft-clipped.
    """
    cdef cystr base
    cdef int pos, readPos
    print("Trying to merge: %s, %s" % (str(pos1), str(pos2)))
    if(pos1.operation != pos2.operation):
        if(pos2.operation == 83):  # if pos2.operation is "S"
            if(pos1.operation == 77):  # if pos1.operation is "M"
                return LayoutPos(pos1.pos, pos1.readPos, pos1.base,
                                 77, pos1.quality, pos1.agreement)
            else:
                return LayoutPos(pos1.pos, pos1.readPos, 78, 78, 0, 0)
                # Returns N for both the cigar op and the base.
                # Quality 0, agreement 0.
        raise ThisIsMadness("Looks like merging these two "
                            "positions just doesn't work (discordance).")
    if(pos1.base == pos2.base):
        print("Agreed base: %s" % chrDict[pos1.base])
        print("Q1: %s. Q2: %s" % (pos1.quality, pos2.quality))
        return LayoutPos(pos1.pos, pos1.readPos, pos1.base, pos1.operation,
                         pos1.quality + pos2.quality,
                         pos1.agreement + pos2.agreement)
    elif(pos1.quality > pos2.quality):
        print("Disagreed! Q1: %s. Q2: %s. " % (pos1.quality, pos2.quality) +
              "B1: %s. B2: %s." % (chr(pos1.base), chr(pos2.base)))
        return LayoutPos(pos1.pos, pos1.readPos, pos1.base, pos1.operation,
                         pos1.quality - pos2.quality, pos1.agreement)
    print("Agreed base: %s" % chrDict[pos1.base])
    print("Q1: %s. Q2: %s" % (pos1.quality, pos2.quality))
    return LayoutPos(pos1.pos, pos1.readPos, pos2.base, pos1.operation,
                     pos2.quality - pos1.quality, pos2.agreement)


cdef tuple MergeLayoutsToList_(Layout_t L1, Layout_t L2):
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
        return [], False
    if(L1.getRefPosForFirstPos() > L2.getRefPosForFirstPos()):
        tmpPos = L1
        L1 = L2
        del tmpPos
    # L1, L2 = sorted((L1, L2), key=omcfp) previous python code
    # Rewritten to avoid the python object omcfp
    offset = L2.getRefPosForFirstPos() - L1.getRefPosForFirstPos()
    print("offset: %s" % offset)
    try:
        return (L1[:offset] + [MergePositions(pos1, pos2) for
                               pos1, pos2 in izip(L1[offset:], L2)] +
                L2[len(L1) - offset:]), True
        '''
        return (L1[:offset] +
                [MergePositions(pos1, pos2) for
                 pos1, pos2 in izip(L1[offset:], L2)] +
                L2[len(L1) - offset:]), True
        '''
    except ThisIsMadness:
        print("ThisIsMadness got thrown!")
        return L1[:] + L2[len(L1) - offset:], False


cpdef Layout_t MergeLayoutsToLayout(Layout_t L1, Layout_t L2):
    cdef list layoutList
    cdef cystr Name
    cdef cython.bint Success
    layoutList, Success = MergeLayoutsToList_(L1, L2)
    if(Success is False):
        L1.tagDict["MP"] = BamTag("MP", tagtype="Z", value="F")
        L2.tagDict["MP"] = BamTag("MP", tagtype="Z", value="F")
        return None
    L1.positions = layoutList
    L1.tagDict["MP"] = BamTag("MP", tagtype="Z", value="T")
    L1.isMerged = True
    L1.update()
    return L1


cpdef cython.bint LayoutsOverlap(Layout_t L1, Layout_t L2):
    return cReadsOverlap(L1.read, L2.read)
