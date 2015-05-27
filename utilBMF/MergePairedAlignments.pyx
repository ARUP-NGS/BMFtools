# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True
from __future__ import division
import pysam
import cython
import numpy as np
from utilBMF.HTSUtils import CigarDict
from operator import attrgetter as oag

cimport pysam.calignmentfile
cimport cython
cimport numpy as np

"""
Contains utilities for merging a pair of overlapping alignments into a single read.
"""
'''
@cython.returns(tuple)
def ProcLayoutArgCigarOp(tuple alignedPair):
    if(alignedPair[1] is None):
'''

@cython.returns(list)
def CigarOpToLayoutPosList(cython.int offset, tuple cigarOp,
                           pysam.calignmentfile.AlignedSegment rec):
    cdef tuple x
    return [LayoutPos(x[1], x[0], CigarDict[cigarOp[0]]) for
            x in rec.aligned_pairs if x[1]]

@cython.returns(tuple)
def makeLayoutTuple(pysam.calignmentfile.AlignedSegment rec):
    cdef tuple pair
    cdef cython.int offset, firstMapped
    cdef list PosList
    offset = 0
    firstMapped = -1
    PosList = []
    for pair in rec.cigar:
        PosList += CigarOpToLayoutPosList(offset, pair, rec)
        if(CigarDict[pair[0]] == "M"):
            firstMapped = offset + rec.pos # Facilitates coordinating merging pairs.
        offset += pair[1]
    return (rec, PosList, firstMapped)


@cython.returns(Layout)
def makeLayout(pysam.calignmentfile.AlignedSegment rec):
    cdef tuple pair
    cdef cython.int offset, firstMapped
    cdef list PosList
    offset = 0
    firstMapped = -1
    PosList = []
    for pair in rec.cigar:
        PosList += CigarOpToLayoutPosList(offset, pair, rec)
        if(CigarDict[pair[0]] == "M"):
            firstMapped = offset  # Facilitates coordinating merging pairs.
        offset += pair[1]
    return Layout(rec, PosList, firstMapped)


@cython.returns(cython.bint)
def lambda1None(tuple i):
    """
    If a tuple in a cigar returns true for this function,
    then that base is either deleted or soft-clipped, which
    you can tell based on whether it is in the middle or the end of a read.
    """
    return i[1] is None


@cython.returns(cython.int)
def getFirstMapped(pysam.calignmentfile.AlignedSegment rec):
    cdef tuple i
    return [i for i in rec.aligned_pairs if i[1] is not None][0][1]
    

cdef class LayoutPos(object):
    """
    Holds one layout position - either part of the reference,
    part of the read, or both.
    """
    def __init__(self, cython.int pos=-1, cython.int readPos=-1,
                 cython.str operation=None):
        self.pos = pos
        self.readPos = readPos
        self.operation = operation


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

oagsk = oag("firstMapped")

def LayoutSortKey(x):
    return oag("firstMapped")(x)


def LayoutSortKey(x):
    return oagsk(x)


@cython.returns(list)
def MergeLayouts(Layout L1, Layout L2):
    L1, L2 = sorted((L1, L2), key=LayoutSortKey)  # First in pair goes first
    pass