# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True
from __future__ import division
import pysam
import cython
import numpy as np
import logging
from utilBMF.HTSUtils import (CigarDict, printlog as pl, PysamToChrDict,
                              ph2chr, TagTypeDict, BamTag)
from operator import attrgetter as oag
from itertools import izip, groupby
from sys import maxint
oagsk = oag("firstMapped")
oagbase = oag("base")
oagop = oag("operation")
oagqual = oag("quality")
oagag = oag("agreement")

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
        self.base = base
        self.quality = quality
        self.agreement = agreement

    def __str__(self):
        return "%s|%s|%s|%s|%s|%s" % (self.pos, self.readPos, self.base,
                                      self.operation, self.quality,
                                      self.agreement)


cdef class Layout(object):
    """
    Holds a read and its layout information.
    """

    @classmethod
    def fromread(cls, pysam.calignmentfile.AlignedSegment rec):
        return cls(*makeLayoutTuple(rec))

    def __getitem__(self, index):
        return self.positions[index]

    def __len__(self):
        return len(self.positions)

    cpdef list get_tags(self):
        return list(self.tagDict.iteritems())

    cpdef cython.str getSeq(self):
        cdef LayoutPos i
        return "".join([i for i in map(oagbase, self.positions) if
                        i not in ["S", "D"]])

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
        return map(ph2chr, self.getQual())

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
        self.tagDict["PV"] = ",".join(map(str, self.getQual()))
        self.tagDict["FA"] = ",".join(map(str, self.getAgreement()))

    def update(self):
        if(self.isMerged):
            # Update it for the merged world!
            # Original template length
            self.tagDict["ot"] = BamTag(("ot", "i", self.tlen))
            # Original mate position
            self.tagDict["op"] = BamTag(("op", "i", self.pnext))
            # Original mapping quality
            self.tagDict["om"] = BamTag(("om", "i", self.mapq))
            self.tlen = 0
            self.pnext = 0
            self.mapq = -1 * maxint
            self.tagDict["MP"] = BamTag(("MP", "A", "T"))
            self.rnext = "*"
        self.update_tags()

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
        return "".join([op + num for op, num in izip(ops, nums)])

    @cython.returns(cython.str)
    def __str__(self):
        """
        Converts the record into a SAM record.
        """
        self.update()
        return "\t".join(map(
                str, [self.Name, self.flag, self.contig,
                      self.getRefPosForFirstPos() + 1, self.mapq,
                      self.getCigarString(), self.rnext, self.pnext,
                      self.tlen, self.getSeq(), self.getQualString()] +
                             map(str, self.get_tags())))

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
        self.tagDict = {tag[0]: BamTag(tag) for tag
                        in rec.get_tags() if tag[0] not in ["PV", "FA"]}
        self.rnext = PysamToChrDict[rec.mrnm]
        self.pnext = rec.mpos
        self.tlen = rec.tlen
        self.isMerged = (rec.has_tag("MP") and rec.opt("MP") == "T")
        self.update()


def LayoutSortKeySK(x, oagsk=oagsk):
    """
    Pre-allocating the sort key function is 20% as fast as allocating
    it at each call. Total sort is ~10% faster by preallocating.

    Faster yet, though, is calling oagsk directly, which is a little over
    an additional 20% speed increase for its use with map.
    """
    return oagsk(x)


@cython.returns(list)
def MergeLayouts(Layout L1, Layout L2, oagsk=oagsk):
    L1, L2 = sorted((L1, L2), key=oagsk)  # First in pair goes first
    pass
