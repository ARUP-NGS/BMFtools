# cython: c_string_type=str, c_string_encoding=ascii
# distutils: sources="include/MPA.c"
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
from utilBMF.ErrorHandling import ThisIsMadness, ImproperArgumentError
from utilBMF.HTSUtils import printlog as pl, BamTag, TrimExt

##DEFINES

oagsk = oag("firstMapped")
oagbase = oag("base")
oagop = oag("operation")
oagqual = oag("quality")
oagag = oag("agreement")
oagtag = oag("tag")
chrDict = {x: chr(x) for x in xrange(33, 126)}


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


cdef class Layout:
    """
    Represents a layout as an array of arrays.
    7 integers each for:
        pos
        readPos
        operation
        base
        quality
        agreement
        mergeAgreed
    , respectively.
    Additional fields:
        firstMapped
            layout index for first "M" base.

    Test for whether a merge has been attempted
    by seeing if mergeAgreed == 1. 1 is unset, 2 is success, 0 is failure.
    """

    def __cinit__(self, AlignedSegment_t read):
        cdef int i, x0, x1
        cdef LayoutLen = 0
        for x0, x1 in read.cigar:
            LayoutLen += x1
        self.Layout = ArrayLayout(layouts=<ArrayLayoutPos_t *>malloc(
            LayoutLen * (sizeof(ArrayLayoutPos_t))), length=LayoutLen)

    def __init__(self, AlignedSegment_t read):
        cdef char CigarOp
        cdef int tmpInt = 0
        cdef int offset = 0
        cdef int CigarOpLen
        cdef tuple tmpTup
        cdef ndarray[int, ndim=1] quals, agrees
        cdef ArrayLayoutPos_t tmpPos
        self.mapq = read.mapq
        self.firstMapped = -1
        try:
            quals = np.array(read.opt("PV").split(","), dtype=np.int32)
        except KeyError:
            pl("Watch out - PV tag not set.", level=logging.DEBUG)
            quals = np.array(read.query_qualities, dtype=np.int32)
            # Let's make sure that these don't need reversal, too!
        try:
            agrees = np.array(read.opt("FA").split(","), dtype=np.int32)
        except KeyError:
            pl("Watch out - FA tag not set.", level=logging.DEBUG)
            agrees = np.array([1] * len(read.sequence), dtype=np.int32)

        # Copy out original alignment information
        self.InitPos = read.pos
        self.tlen = read.tlen
        self.pnext = read.pnext
        self.flag = read.flag
        cdef int readPos = 0
        for tmpTup in read.cigar:
            CigarOp = tmpTup[0]
            CigarOpLen = tmpTup[1]
            for tmpInt in xrange(offset, offset + CigarOpLen):
                if(CigarOp == 0):
                    """
                    Case: 'M'
                    """
                    try:
                        assert read.aligned_pairs[tmpInt][1] is not None
                    except AssertionError:
                        print("Repr of read.align_pairs: %s" % repr(read.aligned_pairs))
                        print("tmpInt: %s." % tmpInt)
                        assert False
                    self.Layout.layouts[tmpInt] = ArrayLayoutPos(
                        pos=read.aligned_pairs[tmpInt][1],
                        readPos=read.aligned_pairs[tmpInt][0],
                        operation=77, base=ord(read.seq[readPos]),
                        quality=quals[readPos], agreement=agrees[readPos],
                        mergeAgreed=1)
                    if(self.firstMapped < 0):
                        self.firstMapped = tmpInt
                    readPos += 1
                elif(CigarOp == 4):
                    """
                    Case: 'S'
                    """
                    try:
                        assert read.aligned_pairs[tmpInt][1] is None
                    except AssertionError:
                        print("Repr of read.align_pairs: %s" % repr(read.aligned_pairs))
                        print("tmpInt: %s." % tmpInt)
                        assert False
                    self.Layout.layouts[tmpInt] = ArrayLayoutPos(
                        pos=-1 * read.pos, readPos=read.aligned_pairs[tmpInt][0],
                        operation=83, base=ord(read.seq[readPos]),
                        quality=quals[readPos], agreement=agrees[readPos],
                        mergeAgreed=1)
                    readPos += 1
                elif(CigarOp == 1):
                    """
                    Case: 'I'
                    """
                    self.Layout.layouts[tmpInt] = ArrayLayoutPos(
                        pos=-1, readPos=read.aligned_pairs[tmpInt][0],
                        operation=73, base=ord(read.seq[readPos]),
                        quality=quals[readPos], agreement=agrees[readPos],
                        mergeAgreed=1)
                    '''
                    self.Layout.layouts[tmpInt].pos = -1
                    self.Layout.layouts[tmpInt].readPos = 
                    self.Layout.layouts[tmpInt].operation = 73
                    self.Layout.layouts[tmpInt].base = ord(read.seq[readPos])
                    self.Layout.layouts[tmpInt].quality = quals[readPos]
                    self.Layout.layouts[tmpInt].agreement = agrees[readPos]
                    self.Layout.layouts[tmpInt].mergeAgreed = 1
                    '''
                    readPos += 1
                elif(CigarOp == 2):
                    """
                    Case: 'D'
                    """
                    self.Layout.layouts[tmpInt] = ArrayLayoutPos(
                        pos=read.aligned_pairs[tmpInt][1],
                        readPos=-1,
                        operation=68, base=68,
                        quality=-1, agreement=-1, mergeAgreed=1)
                else:
                    raise NotImplementedError(
                        "Only MIDS cigar operations currently supported. If "
                        "you have an application that could use further "
                        "support, please contact me.")
            offset += CigarOpLen

    cdef bint cPosIsMapped(self, int position):
        return self.Layout.layouts[position].operation == 77  # == "M"

    cpdef bint posIsMapped(self, int position):
        return self.cPosIsMapped(position)

    cdef int getFirstMappedReadPos(self):
        cdef int i
        for i in range(self.Layout.length):
            if(self.Layout.layouts[i].operation == 77):
                return i

    cdef int getFirstAlignedRefPos(self):
        return getFirstAlignedRefPos(self.Layout)

    cdef int getFirstMappedRefPos(self):
        cdef int tmpInt
        for tmpInt in range(self.Layout.length):
            if(self.Layout.layouts[tmpInt].operation == 77):
                return self.Layout.layouts[tmpInt].pos
            # Operation is M, returns the ref position.
        raise ImproperArgumentError(
            "ArrayLayout has no 'M' cigar operation positions. "
            "This read can't be layed out???")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef ndarray[np.int16_t, ndim=1] cGetQual(self):
        # Ask if >= 0. My tests say it's ~1% faster to ask (> -1) than (>= 0).
        # pos.base == 66 for "B", which is a blank spot.
        cdef int tmpInt
        return np.array([self.Layout.layouts[tmpInt].quality for
                         tmpInt in range(self.Layout.length) if
                         self.Layout.layouts[tmpInt].operation != 68],
                        dtype=np.int16)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef ndarray[int, ndim=1] getQual(self):
        return self.cGetQual()

    cdef cystr cGetQualString(self):
        cdef int i
        return "".join([ph2chrDict[i] for i in self.cGetQual()])

    cpdef cystr getQualString(self):
        return self.cGetQualString()

    def __dealloc__(self):
        free(self.Layout.layouts)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef ndarray[char] cGetSeqArr(self):
        """Returns a character array of the base calls
        if the operations aren't "D" (68)
        """
        cdef int i
        return np.char.array([chrDict[self.Layout.layouts[i].base]
                              for i in xrange(self.Layout.length) if
                              self.Layout.layouts[i].operation != 68],
                             itemsize=1)

    cpdef ndarray[char] getSeqArr(self):
        return self.cGetSeqArr()

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef cystr getSeq(self):
        cdef int i
        return "".join([chrDict[self.Layout.layouts[i].base] for i in
                        xrange(self.Layout.length) if
                        self.Layout.layouts[i].operation != 68])

    @cython.returns(cystr)
    def __str__(self):
        cdef int i
        return "#".join(map(ALPToStr, [self.Layout.layouts[i] for
                                       i in range(self.Layout.length)]))

    cdef MergeLayouts_in_place(self, ArrayLayout_t pairedLayout):
        self.Layout = MergeLayouts(self.Layout, pairedLayout)

'''

Dinosaur code from before the struct version of Layout. Keeping around since
I'll have to recode most of this for the new object scheme.
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

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef ndarray[char] getSeqArr(self):
        """Returns a character array of the base calls
        if the base calls aren't "D" (68)
        """
        cdef char i
        cdef LayoutPos_t pos
        return np.char.array([chrDict[pos.base]
                              for pos in self.positions if
                              pos.operation != 68 and
                              pos.agreement > -1],
                             itemsize=1)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef cystr getSeq(self):
        return self.getSeqArr().tostring()

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int cGetRefPosForFirstPos(self):
        """cdef class wrapped by getRefPosForFirstPos
        """
        cdef LayoutPos_t i
        cdef int count
        for count, i in enumerate(self):
            if(i.operation == "M"):
                return i.pos - count

    cpdef int getRefPosForFirstPos(self):
        return self.cGetRefPosForFirstPos()

    @cython.boundscheck(False)
    @cython.wraparound(False)
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
                        dtype=np.int32)

    cpdef ndarray[int, ndim=1] getQual(self):
        return self.cGetQual()

    cdef cystr cGetQualString(self):
        cdef int i
        return "".join([ph2chrDict[i] for i in self.cGetQual()])

    cpdef cystr getQualString(self):
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

    cpdef ndarray[np.int16_t, ndim=1] getMergedPositions(self):
        """
        Returns the indices of the bases in the merged read
        which have been merged.
        """
        return self.cGetMergedPositions()

    cpdef list getPositions(self):
        cdef LayoutPos_t tmpPos
        return [tmpPos for tmpPos in self.positions]

    cpdef ndarray[char, ndim=1] getMergeAgreements(self):
        """
        Returns the indices of the bases in the merged read
        which have been merged successfully.
        """
        return self.cGetMergeAgreements()

    cdef ndarray[np.int16_t, ndim=1] cGetMergedPositions(self):
        """
        Returns the indices of the bases in the merged read
        which have been merged.
        """
        cdef LayoutPos_t lpos
        return np.array([lpos.readPos for lpos in self.positions if
                         lpos.merged],
                        dtype=np.int16)

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
            self.tagDict["MP"] = BamTag("MP", "Z", "T")

    cpdef update_tags(self):
        self.update_tags_()

    def update(self):
        cdef LayoutPos_t pos
        cdef int count
        self.update_tags()
        if(self.isMerged):
            self.tlen = len(self.getSeqArr())
            self.pnext = 0
            # Only change the original mapq to -1 if the tagDict entry om is
            # present. (Original Mapping)
            try:
                self.tagDict["om"]
                pass
            except KeyError:
                self.mapq = -1
            self.rnext = "*"
            self.flag = 2 + (16 if(self.is_reverse) else 32)
            for count, pos in enumerate(self):
                pos.readPos = count

    @cython.returns(ndarray)
    def getOperations(self, oagop=oagop):
        """
        In [15]: %timeit [chr(p.operation) for p in l1.positions]
        100000 loops, best of 3: 9.35 us per loop

        In [16]: %timeit map(chr, [p.operation for p in l1.positions])
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

    @cython.returns(AlignedSegment_t)
    def __read__(self):
        cdef AlignedSegment_t newRead
        cdef int i
        cdef BamTag_t BT
        newRead = pysam.AlignedSegment()
        newRead.query_name = self.Name
        newRead.flag = self.getFlag()
        newRead.reference_id = PysamToChrDict[self.contig]
        newRead.query_sequence = self.getSeq()
        print(str(self.read))
        newRead.query_qualities = [93 if(i > 92) else i for
                                   i in self.cGetQual()]
        newRead.reference_start = self.getAlignmentStart()
        newRead.cigarstring = self.getCigarString()
        newRead.tlen = self.tlen
        newRead.mapq = self.mapq
        newRead.tags = [(BT.tag, BT.value) for BT in self.tagDict.itervalues()]
        newRead.next_reference_id = PysamToChrDict[self.rnext]
        newRead.pnext = self.pnext
        return newRead

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
'''


def MergePairedAlignments(cystr inBAM, cystr outBAM=None,
                          bint pipe=False, int readLength=-1):
    cdef AlignedSegment_t read, read1, read2
    cdef Layout_t Layout1, Layout2, retLayout
    cdef int rLen2
    cdef int count = 0
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    if(readLength < 0):
        pl("readLength not set - inferring.")
        readLength = len(inHandle.next().seq)
        inHandle = pysam.AlignmentFile(inBAM, "rb")
    rLen2 = 2 * readLength  # If tlen >= rLen2, no overlap.
    if(outBAM is None):
        outBAM = TrimExt(inBAM) + ".PairMergeProcessed.bam"
    if(pipe):
        outHandle = pysam.AlignmentFile("-", "wb",
                                        header=inHandle.header)
    else:
        outHandle = pysam.AlignmentFile(outBAM, "wb",
                                        header=inHandle.header)
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
        if(read1.reference_id != read2.reference_id or read1.is_unmapped or
           read2.is_unmapped or abs(read1.tlen) >= rLen2):
            ohw(read1)
            ohw(read2)
            continue
        try:
            assert(read1.query_name == read2.query_name)
        except AssertionError:
            raise ThisIsMadness("Bam is either not name sorted or you are "
                                "missing a read from the pair around read "
                                "# %s in the bam." % count)
        '''
        Layout1 = Layout(read1)
        Layout2 = Layout(read2)
        retLayout = MergeLayoutsToLayout(Layout1, Layout2)
        if(retLayout is None):
            read1.setTag("MP", "F")
            read2.setTag("MP", "F")
            ohw(read1)
            ohw(read2)
            continue
        ohw(retLayout.__read__())
        '''
    outHandle.close()
    inHandle.close()
    return outBAM


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int getLayoutLen(AlignedSegment_t read):
    cdef tuple tmpTup
    cdef int lensum = 0
    if(read.cigarstring is None):
        raise ImproperArgumentError(
            "read " + read.query_name +
            " is unmapped - no such thing as a layout length!")
    for tmpTup in read.cigar:
        lensum += tmpTup[1]
    return lensum


cpdef cystr ALPToStr(ArrayLayoutPos_t ALP):
    return "|".join(map(str, [ALP.pos, ALP.readPos, ALP.quality,
                              ALP.agreement, chr(ALP.operation),
                              chr(ALP.base), ALP.mergeAgreed]))
