# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True
from __future__ import division
import pysam
import cython
import numpy as np
import logging
from utilBMF.HTSUtils import (CigarDict, printlog as pl, PysamToChrDict,
                              ph2chr, TagTypeDict, BamTag)
from utilBMF.ErrorHandling import ThisIsMadness
from operator import attrgetter as oag, methodcaller as mc
from itertools import izip, groupby
from sys import maxint
oagsk = oag("firstMapped")
omcfp = mc("getRefPosForFirstPos")
oagbase = oag("base")
oagop = oag("operation")
oagqual = oag("quality")
oagag = oag("agreement")
oagtag = oag("tag")

cimport pysam.calignmentfile
cimport cython
cimport numpy as np

"""
Contains utilities for merging a pair of overlapping alignments
into a single read.
"""


@cython.returns(list)
def CigarOpToLayoutPosList(cython.int offset, tuple cigarOp,
                           pysam.calignmentfile.AlignedSegment rec):
    cdef tuple x
    cdef cython.str CigarChar
    cdef np.ndarray[cython.long, ndim = 1] quals, agrees
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
        quals = np.array(rec.quality, dtype=np.int64)
    try:
        agrees = np.array(rec.opt("FA").split(","), dtype=np.int64)
    except KeyError:
        pl("Watch out - FA tag not set.", level=logging.DEBUG)
        agrees = np.array([1] * len(rec.sequence), dtype=np.int64)
    return [LayoutPos(pos=x[1], readPos=x[0], operation=CigarChar,
                      base=rec.seq[x[0]], quality=quals[x[0]],
                      agreement=agrees[x[0]]) if
            x[1] is not None and x[0] is not None else
            LayoutPos(-1, x[0], operation=CigarChar, base=rec.seq[x[0]],
                      quality=quals[x[0]], agreement=agrees[x[0]]) if
            x[0] is not None else
            LayoutPos(x[1], -1, operation=CigarChar, base=CigarChar,
                      quality=quals[x[0]], agreement=agrees[x[0]]) if
            CigarChar == "S" else
            LayoutPos(x[1], -1, operation=CigarChar, base=CigarChar,
                      quality=-1, agreement=-1)
            for x in rec.aligned_pairs[offset:offset + cigarOp[1]]]


@cython.returns(Layout)
def makeLayout(pysam.calignmentfile.AlignedSegment rec):
    return Layout(makeLayoutTuple)


@cython.returns(tuple)
def makeLayoutTuple(pysam.calignmentfile.AlignedSegment rec):
    cdef tuple pair
    cdef cython.int offset, firstMapped, i
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


@cython.returns(cython.int)
def getFirstMappedRefPos(pysam.calignmentfile.AlignedSegment rec):
    cdef tuple i
    return [i for i in rec.aligned_pairs if i[1] is not None][0][1]


cdef class LayoutPos(object):
    """
    Holds one layout position - either part of the reference,
    part of the read, or both.

    """
    def __init__(self, cython.int pos=-1, cython.int readPos=-1,
                 cython.str base=None, cython.str operation=None,
                 cython.int quality=-1, cython.int agreement=-1):
        self.pos = pos
        self.readPos = readPos
        self.operation = operation
        self.base = base if(quality > 2) else "N"
        self.quality = quality if(base != "N") else 0
        self.agreement = agreement

    def __str__(self):
        return "%s|%s|%s|%s|%s|%s" % (self.pos, self.readPos, self.base,
                                      self.operation, self.quality,
                                      self.agreement)


cdef class Layout(object):
    """
    Holds a read and its layout information.

    This doctest was written so that it would load in one read,
    make the string, and then hash that value. Since it wouldn't
    match the interpreter's output to have a gigantic line and it would have
    to violate pep8, I decided to test the value by its hash rather than by
    string agreement.
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

    cpdef cython.str getSeq(self):
        cdef cython.str i
        return "".join([i for i in map(oagbase, self.positions) if
                        i != "S" and i != "D"])

    @cython.returns(cython.int)
    def getRefPosForFirstPos(self):
        cdef LayoutPos i
        cdef cython.int count
        for count, i in enumerate(self):
            if(i.operation == "M"):
                return i.pos - count

    @cython.returns(list)
    def getAgreement(self, oagag=oagag):
        cdef cython.int i
        return [i for i in map(oagag, self.positions) if i >= 0]

    @cython.returns(list)
    def getQual(self, oagqual=oagqual):
        cdef cython.int i
        return [i for i in map(oagqual, self.positions) if i >= 0]

    def getQualString(self):
        return "".join(map(ph2chr, self.getQual()))

    @cython.returns(cython.int)
    def getLastRefPos(self):
        """
        Finds the reference position which "matches" the last base in the
        layout. This counts for any cigar operation.
        """
        cdef LayoutPos i
        cdef cython.int lastMPos, countFromMPos
        lastMPos = 0
        for count, i in enumerate(self):
            countFromMPos += 1
            if(i.operation == "M"):
                lastMPos = i.pos
                countFromMPos = 0
        return lastMPos + countFromMPos

    def update_tags(self):
        self.tagDict["PV"] = BamTag("PV", "Z",
                                    ",".join(map(str, self.getQual())))
        self.tagDict["FA"] = BamTag("FA", "Z",
                                    ",".join(map(str, self.getAgreement())))

    def update_flag(self):
        """
        TODO: Fill this in!
        """
        pass

    def update(self):
        cdef LayoutPos pos
        cdef cython.int count
        if(self.isMerged):
            # Update it for the merged world!
            # Original template length
            self.tagDict["ot"] = BamTag("ot", "i", self.tlen)
            # Original mate position
            self.tagDict["mp"] = BamTag("mp", "i", self.pnext)
            # Original mapping quality
            self.tagDict["om"] = BamTag("om", "i", self.mapq)
            self.tagDict["op"] = BamTag("op", "i", self.InitPos)
            self.tlen = 0
            self.pnext = 0
            self.mapq = -1 * maxint
            self.tagDict["MP"] = BamTag("MP", "A", "T")
            self.rnext = "*"
            for count, pos in enumerate(self):
                pos.readPos = count
        self.update_tags()
        self.update_flag()

    @cython.returns(list)
    def getOperations(self, oagop=oagop):
        return map(oagop, self.positions)

    cpdef cython.str getCigarString(self):
        cdef cython.str op, num
        cdef list nums = []
        cdef list ops = []
        for k, g in groupby(self.getOperations()):
            nums.append(str(len(list(g))))
            ops.append(k)
        return "".join([num + op for op, num in izip(ops, nums)])

    @cython.returns(list)
    def get_tags(self, oagtag=oagtag):
        self.update_tags()
        return sorted(self.tagDict.itervalues(), key=oagtag)

    @cython.returns(cython.str)
    def __str__(self):
        """
        Converts the record into a SAM record.
        Note: the position is incremented by 1 because SAM positions are
        1-based instead of 0-based.
        """
        self.update()
        return "\t".join(map(
                str, [self.Name, self.flag, self.contig,
                      self.getRefPosForFirstPos() + 1, self.mapq,
                      self.getCigarString(), self.rnext, self.pnext + 1,
                      self.tlen, self.getSeq(), self.getQualString()] +
                             self.get_tags()))

    def __init__(self, pysam.calignmentfile.AlignedSegment rec,
                 list layoutPositions, cython.int firstMapped,
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
        self.tagDict = {tag[0]: BamTag.fromtuple(tag) for tag
                        in rec.get_tags() if tag[0] not in ["PV", "FA"]}
        self.rnext = PysamToChrDict[rec.mrnm]
        self.pnext = rec.mpos
        self.tlen = rec.tlen
        self.isMerged = (rec.has_tag("MP") and rec.opt("MP") == "T")


def LayoutSortKeySK(x, oagsk=oagsk):
    """
    Pre-allocating the sort key function is 20% as fast as allocating
    it at each call. Total sort is ~10% faster by preallocating.

    Faster yet, though, is calling oagsk directly, which is a little over
    an additional 20% speed increase for its use with map.
    """
    return oagsk(x)

@cython.returns(LayoutPos)
def MergePositions(LayoutPos pos1, LayoutPos pos2):
    cdef cython.str base, operation
    cdef cython.int pos, readPos
    if(pos1.operation != pos2.operation):
        raise ThisIsMadness("Looks like merging these two "
                            "positions just doesn't work (discordance).")
    if(pos1.base == pos2.base):
        return LayoutPos(pos1.pos, pos1.readPos, pos1.base, pos1.operation,
                         pos1.quality + pos2.quality,
                         pos1.agreement + pos2.agreement)
    elif(pos1.quality > pos2.quality):
        return LayoutPos(pos1.pos, pos1.readPos, pos1.base, pos1.operation,
                         pos1.quality - pos2.quality, pos1.agreement)
    return LayoutPos(pos1.pos, pos1.readPos, pos2.base, pos1.operation,
                     pos2.quality - pos1.quality, pos2.agreement)

@cython.returns(tuple)
def MergeLayoutsToList(Layout L1, Layout L2, omcfp=omcfp):
    """
    Merges two Layouts into a list of layout positions.

    First, it takes the positions before the overlap starts.
    Then it merges the overlap position by position.
    Then it takes the positions after the overlap ends.
    :param Layout L1: One Layout object
    :param Layout L2: Another Layout object
    :param oagsk: A method caller function. Locally redefined for speed.
    
    :return list Merged Positions
    :return bool Whether the merge was successful
    """
    cdef cython.int offset
    L1, L2 = sorted((L1, L2), key=omcfp)  # First in pair goes first
    offset = L2.getRefPosForFirstPos() - L1.getRefPosForFirstPos()
    try:
        return (L1[:offset] +
                [MergePositions(pos1, pos2) for
                 pos1, pos2 in izip(L1[offset:], L2)] +
                L2[len(L1) - offset:]), True
    except ThisIsMadness:
        return L1[:] + L2[len(L1) - offset:], False
