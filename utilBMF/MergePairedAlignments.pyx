# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True
from __future__ import division
import pysam
import cython
import numpy as np
import logging
from utilBMF.HTSUtils import CigarDict, printlog as pl
from operator import attrgetter as oag
oagsk = oag("firstMapped")
oagbase = oag("base")

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
    def __init__(self, pysam.calignmentfile.AlignedSegment rec,
                 list layoutPositions, cython.int firstMapped):
        self.read = rec
        self.positions = layoutPositions
        self.firstMapped = firstMapped

    @classmethod
    def fromread(cls, pysam.calignmentfile.AlignedSegment rec):
        return cls(*makeLayoutTuple(rec))

    def __getitem__(self, index):
        return self.positions[index]

    @cython.returns(cython.str)
    def getSeq(self):
        cdef LayoutPos i
        return "".join([i for i in map(oagbase, self.positions) if
                        i not in ["S", "D"]])


def LayoutSortKeySK(x):
    """
    Pre-allocating the sort key function is 20% as fast as allocating
    it at each call. Total sort is ~10% faster by preallocating.

    Faster yet, though, is calling oagsk directly, which is a little over
    an additional 20% speed increase for its use with map.
    """
    return oagsk(x)


@cython.returns(list)
def MergeLayouts(Layout L1, Layout L2):
    L1, L2 = sorted((L1, L2), key=oagsk)  # First in pair goes first
    pass