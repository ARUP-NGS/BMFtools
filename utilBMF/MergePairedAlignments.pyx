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
from pysam import AlignedSegment

## Standard Library local imports
from array import array
from itertools import izip, groupby
from operator import attrgetter as oag, methodcaller as mc, itemgetter as oig
from subprocess import check_call
from sys import maxint
try:
    from re2 import split as rsplit
except ImportError:
    from re import split as rsplit

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
oig0 = oig(0)
oig1 = oig(1)
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
CigarStrDict = {77: "M", 73: "I", 68: "D", 78: "N", 83: "S", 72: "H",
                80: "P", 61: "=", 88: "X"}

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
def getFirstMappedRefPos(AlignedSegment_t rec):
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
        # C Definitions
        cdef char CigarOp, CigarOpLen
        cdef int tmpInt = 0
        cdef int offset = 0
        cdef int readPos = 0
        cdef ndarray[int, ndim=1] quals, agrees
        cdef ArrayLayoutPos_t tmpPos
        cdef tuple tag

        # Parsing out base-by-base information.
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
        self.pos = read.pos
        self.tlen = read.tlen
        self.pnext = read.pnext
        self.flag = read.flag
        self.Name = read.query_name
        self.reference_id = read.reference_id
        self.is_reverse = read.is_reverse
        self.tagDict = {tag[0]: BamTag.fromtuple(tag) for tag
                        in read.get_tags() if tag[0] not in ["PV", "FA"]}
        self.isMerged = (read.has_tag("MP") and read.opt("MP") == "T")

        # Parse out the read cigar.
        for CigarOp, CigarOpLen in read.cigar:
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
        cdef size_t i
        for i in range(self.Layout.length):
            if(self.Layout.layouts[i].operation == 77):
                return i

    cdef int getLastAlignedRefPos(self):
        cdef size_t i, counter
        counter = 0
        for i from self.Layout.length > i >= 0:
            if(self.Layout.layouts[i].operation == 77):
                return self.Layout.layouts[i].pos + counter
            counter += 1

    cdef int getFirstAlignedRefPos(self):
        return getFirstAlignedRefPos(self.Layout)

    cdef int getAlignmentStart(self):
        cdef int tmpInt
        for tmpInt in range(self.Layout.length):
            if(self.Layout.layouts[tmpInt].operation == 77):
                return self.Layout.layouts[tmpInt].pos
            # Operation is M, returns the ref position.
        raise ImproperArgumentError(
            "ArrayLayout has no 'M' cigar operation positions. "
            "This read can't be layed out???")

    cdef ndarray[int, ndim=1] cGetQual(self):
        # Ask if >= 0. My tests say it's ~1% faster to ask (> -1) than (>= 0).
        # pos.base == 66 for "B", which is a blank spot.
        cdef int tmpInt
        return np.array([self.Layout.layouts[tmpInt].quality for
                         tmpInt in range(self.Layout.length) if
                         self.Layout.layouts[tmpInt].operation != 68],
                        dtype=np.int32)

    cdef ndarray[int, ndim=1] cGetAgreementSlice(self, size_t start=0):
        return np.array([self.Layout.layouts[tmpInt].agreement for tmpInt
                         in range(start, self.Layout.length) if
                         self.Layout.layouts[tmpInt].operation != 68])

    cdef ndarray[int, ndim=1] cGetQualSlice(self, size_t start=0):
        return np.array([self.Layout.layouts[tmpInt].quality for tmpInt
                         in range(start, self.Layout.length) if
                         self.Layout.layouts[tmpInt].operation != 68])

    cpdef ndarray[int, ndim=1] getQual(self):
        return self.cGetQual()

    cdef cystr cGetQualStringSlice(self, size_t start=0):
        cdef int i
        return "".join([ph2chrDict[i] for i
                        in self.cGetQualSlice(start)])

    cdef cystr cGetQualString(self):
        cdef int i
        return "".join([ph2chrDict[i] for i in self.cGetQual()])

    cpdef cystr getQualString(self):
        return self.cGetQualString()

    def __dealloc__(self):
        free(self.Layout.layouts)

    cdef ndarray[char] cGetSeqArrSlice(self, size_t start=0):
        """Returns a character array of the base calls
        if the operations aren't "D" (68)
        """
        cdef size_t i
        return np.char.array([chrDict[self.Layout.layouts[i].base]
                              for i in xrange(start, self.Layout.length) if
                              self.Layout.layouts[i].operation != 68],
                             itemsize=1)

    cdef ndarray[char] cGetSeqArr(self):
        """Returns a character array of the base calls
        if the operations aren't "D" (68)
        """
        cdef size_t i
        return np.char.array([chrDict[self.Layout.layouts[i].base]
                              for i in xrange(self.Layout.length) if
                              self.Layout.layouts[i].operation != 68],
                             itemsize=1)

    cpdef ndarray[char] getSeqArr(self):
        return self.cGetSeqArr()

    cpdef cystr getSeq(self):
        cdef int i
        return "".join([chrDict[self.Layout.layouts[i].base] for i in
                        xrange(self.Layout.length) if
                        self.Layout.layouts[i].operation != 68])

    @cython.returns(cystr)
    def __str__(self):
        cdef ArrayLayoutPos_t tmpALP
        return "#".join(map(ALPToStr, [ALP for ALP in
                                       self.Layout.layouts[:self.Layout.length]
                                       ]))
    '''
    cpdef AlignedSegment_t MergeLayoutToAS(self, Layout_t pairedLayout,
                                           AlignedSegment_t template):
        cdef size_t tmpInt
        cdef cystr tmpQual
        cdef char offset
        # cdef AlignedSegment_t retAS = self.read
        retAS = CopyAlignedSegment(template)
        offset = self.MergeLayouts_in_place(pairedLayout.Layout)
        if(self.MergeSuccess is False):
            return None
        retAS.cigarstring += pairedLayout.cGetCigarStringSlice(offset)
        retAS.cigar = FlattenCigar(retAS.cigar)
        tmpQual = retAS.qual
        retAS.seq += pairedLayout.cGetSeqArrSlice(offset).tostring()
        retAS.qual = tmpQual + pairedLayout.cGetQualStringSlice(offset)
        retAS.set_tag("PV", retAS.opt("PV") + "," +
                          ",".join(self.cGetQualSlice(offset).astype(str)))
        retAS.set_tag("FA",
                      (retAS.opt("FA") +
                       "," +
                       ",".join(self.cGetAgreementSlice(offset).astype(str))))
        retAS.set_tags(self.get_merge_tags())
        retAS.mapq = 255  # Meaning mapping quality is not meaningful without a realignment.
        retAS.reference_id = self.reference_id
        retAS.pos = self.getAlignmentStart()
        return retAS
    '''

    cdef char MergeLayouts_in_place(self, ArrayLayout_t pairedLayout):
        print("Begging MergeLayouts_in_place")
        cdef char offset
        offset = MergeOverlappedLayouts(self.Layout, pairedLayout)
        cdef int tmpInt
        self.MergeSuccess = True
        for tmpInt in range(self.Layout.length):
            print("ALP String: %s" % ALPToStr(self.Layout.layouts[tmpInt]))
            if(self.Layout.layouts[tmpInt].mergeAgreed < 0):
                self.MergeSuccess = False
        self.update()
        self.isMerged = True
        return offset

    cdef ndarray[np.int16_t, ndim=1] getMergeAgreements(self):
        cdef int tmpInt
        return np.array([self.Layout.layouts[tmpInt].readPos for tmpInt
                         in range(self.Layout.length) if
                         self.Layout.layouts[tmpInt].mergeAgreed == 2],
                        dtype=np.int16)

    cdef ndarray[np.int16_t, ndim=1] getMergedPositions(self):
        cdef int tmpInt
        return np.array([self.Layout.layouts[tmpInt].readPos for tmpInt
                         in range(self.Layout.length) if
                         self.Layout.layouts[tmpInt].mergeAgreed != 1],
                        dtype=np.int16)

    cdef set_merge_tags_BT(self):
        cdef cystr tag
        cdef object value
        for tag, value in self.get_merge_tags():
            self.tagDict[tag] = BamTag.fromtuple((tag, value))

    cdef list get_merge_tags(self):
        cdef ndarray[np.int32_t, ndim=1] GenDiscPos
        cdef ndarray[np.int16_t, ndim=1] ReadDiscPos
        cdef str PMStr, MAStr, DGStr, DRStr

        assert self.isMerged

        GenDiscPos = self.cGetGenomicDiscordantPositions()
        ReadDiscPos = self.cGetReadDiscordantPositions()
        PMStr = ",".join(self.getMergedPositions().astype(str))
        MAStr = ",".join(self.getMergeAgreements().astype(str))
        DGStr = ",".join(GenDiscPos.astype(str))
        DRStr = ",".join(ReadDiscPos.astype(str))
        if(len(GenDiscPos) > 0):
            return [("PM", PMStr), ("MA", MAStr),
                    ("ot", self.tlen), ("mp", self.pnext),
                    ("om", self.mapq), ("MP", "T"),
                    ("DG", DGStr),("DR", DRStr)]
        return [("PM", PMStr), ("MA", MAStr),
                ("ot", self.tlen), ("mp", self.pnext),
                ("om", self.mapq), ("MP", "T")]

    cpdef cystr __sam__(self, Layout_t pairedLayout):
        return self.__csam__(pairedLayout)

    cdef cystr __csam__(self, Layout_t pairedLayout):
        """
        Converts the record into a SAM record.
        Note: the position is incremented by 1 because SAM positions are
        1-based instead of 0-based.
        """
        cdef cystr NewCigarString
        cdef int offset
        offset = self.MergeLayouts_in_place(pairedLayout.Layout)
        if(self.MergeSuccess is False):
            return None
        self.update()
        self.tagDict["PV"].value += (
            "," + ",".join(self.cGetQualSlice(offset).astype(str)))
        self.tagDict["FA"].value += (
            "," + ",".join(self.cGetAgreementSlice(offset).astype(str)))
        self.tagDict["FM"].value *= 2
        # Double the number of "family members" to describe merging.
        NewCigarString = FlattenCigarString(
            self.getCigarString + pairedLayout.cGetCigarStringSlice(offset))
        '''
        return "\t".join(map(
                str, [self.Name, self.getFlag(), self.contig,
                      self.getAlignmentStart() + 1, self.mapq,
                      self.getCigarString(), self.rnext, self.pnext + 1,
                      self.tlen, self.getSeq(), self.getQualString()] +
                self.get_tag_string()))
        '''
        return (self.Name + "\t%s" % self.getFlag() + "\t" + self.contig +
                "\t%s\t%s" % (self.getAlignmentStart() + 1, self.mapq) +
                "\t" + NewCigarString +
                "\t%s\t%s\t%s\t" % (self.rnext, self.pnext + 1, self.tlen) +
                self.getSeq() +
                pairedLayout.cGetSeqArrSlice(offset).tostring() + "\t" +
                self.getQualString() +
                pairedLayout.cGetQualStringSlice(offset) + "\t" +
                + self.get_tag_string())

    cdef get_tag_string(self):
        cdef BamTag_t BT
        return "\t".join([BT.__str__() for BT in self.tagDict.itervalues()])

    cdef update_read_positions(self):
        cdef ArrayLayoutPos_t tmpALP
        cdef int tmpInt
        for tmpALP in self.Layout.layouts[:self.Layout.length]:
            self.tmpALP.readPos = tmpInt if(tmpALP.operation != 68) else -1
            tmpInt += 1

    def update(self):
        cdef int count
        if(self.isMerged):
            self.set_merge_tags_BT()
            self.tlen = len(self.getSeqArr())
            self.pnext = 0
            self.rnext = self.reference_id
            # Only change the original mapq to -1 if the tagDict entry om is
            # present. (Original Mapping)
            try:
                self.tagDict["om"]
                self.mapq = 255
            except KeyError:
                pass
            self.flag = 2 + (16 if(self.is_reverse) else 32)
            self.update_read_positions()
        # self.update_read()

    '''
    cdef update_read(self):
        cdef int tmpInt
        cdef BamTag_t BT
        cdef list tmpList
        self.read.seq = self.getSeq()
        self.read.query_qualities = [93 if(tmpInt > 92) else 0 if(tmpInt < 0)
                                     else tmpInt for
                                     tmpInt in self.cGetQual()]
        self.read.flag = self.flag
        self.read.tlen = self.tlen
        self.read.pnext = self.pnext
        self.read.pos = self.getAlignmentStart()
        self.read.cigarstring = self.cGetCigarString()
        self.read.mapq = self.mapq
        self.read.tags = sorted([(BT.tag, BT.value) for BT in
                                 self.tagDict.itervalues()],
                                key=oig0)
    '''

    cdef ndarray[int, ndim=1] cGetGenomicDiscordantPositions(self):
        cdef ArrayLayoutPos_t tmpLayoutPos
        return np.array([tmpLayoutPos.pos for tmpLayoutPos in
                         self.Layout.layouts[:self.Layout.length] if
                         tmpLayoutPos.mergeAgreed == 0 and
                         tmpLayoutPos.operation != 66],
                        dtype=np.int32)

    cdef ndarray[np.int16_t, ndim=1] cGetReadDiscordantPositions(self):
        cdef ArrayLayoutPos_t tmpLayoutPos
        return np.array([tmpLayoutPos.readPos for tmpLayoutPos in
                         self.Layout.layouts[:self.Layout.length] if
                         tmpLayoutPos.mergeAgreed == 0 and
                         tmpLayoutPos.operation != 66],
                        dtype=np.int16)

    cpdef ndarray[int, ndim=1] getAgreement(self):
        """cpdef wrapper of cGetAgreement
        """
        return self.cGetAgreement()

    cdef ndarray[uint16_t, ndim=1] cGetAgreement(self):
        cdef ArrayLayoutPos_t tmpLayoutPos
        # Ask if >= 0. My tests say it's ~1% faster to ask (> -1) than (>= 0).
        return np.array([tmpLayoutPos.agreement for tmpLayoutPos in
                         self.Layout.layouts[:self.Layout.length] if
                         tmpLayoutPos.operation != 66],
                        dtype=np.int16)

    cpdef ndarray[np.int8_t, ndim=1] getOperations(self):
        return self.cGetOperations()

    cdef ndarray[np.int8_t, ndim=1] cGetOperationsSlice(self, size_t start=0):
        cdef int tmpInt
        return np.array([self.Layout.layouts[tmpInt].operation for tmpInt in
                         range(start, self.Layout.length)], dtype=np.int8)

    cdef ndarray[np.int8_t, ndim=1] cGetOperations(self):
        cdef int tmpInt
        return np.array([self.Layout.layouts[tmpInt].operation for tmpInt in
                         range(self.Layout.length)], dtype=np.int8)

    cdef cystr cGetCigarStringSlice(self, size_t start=0):
        return "".join([str(len(list(g))) + CigarStrDict[k] for
                        k, g in groupby(self.cGetOperationsSlice(start))])

    cdef cystr cGetCigarString(self):
        return "".join([str(len(list(g))) + CigarStrDict[k] for
                        k, g in groupby(self.cGetOperations())])

    cpdef cystr getCigarString(self):
        return self.cGetCigarString()

    def getFlag(self):
        self.update()
        return self.flag


def MergePairedAlignments(cystr inBAM, cystr outBAM=None,
                          int readLength=-1,
                          cystr outMerge=None):
    cdef AlignedSegment_t read, read1, read2, newRead
    cdef Layout_t Layout1, Layout2, retLayout
    cdef int count = 0
    cdef cystr newSamRead
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    if(readLength < 0):
        pl("readLength not set - inferring.")
        readLength = len(inHandle.next().seq)
        inHandle = pysam.AlignmentFile(inBAM, "rb")
    outHandle = pysam.AlignmentFile("-", "w",
                                    header=inHandle.header)
    outHandleMerge = sys.stdout
    outHandleMerge.write(inHandle.text)
    for read in inHandle:
        print("Name: %s" % read.query_name)
        count += 1
        if(read.is_supplementary or read.is_secondary):
            outHandle.write(read)
            print("I am now continuing over this record in MPA for 2nd/Supp.")
            continue
        if(read.is_read1):
            read1 = read
            print("I am now continuing over this record in MPA because it's read1.")
            continue
        elif(read.is_read2):
            read2 = read
        else:
            print("Read is neither read 1 nor read 2: %s" % (not read.is_read1 and not read.is_read2))
            raise ThisIsMadness("Read is neither read 1 nor read 2? %s" % read.query_name)
        try:
            assert(read1.query_name == read2.query_name)
        except AssertionError:
            print("read1.query_name: %s. read2.query_name: %s" % (read1.query_name, read2.query_name))
            raise ThisIsMadness("Bam is either not name sorted or you are "
                                "missing a read from the pair around read "
                                "# %s in the bam." % count)
        if(cReadsOverlap(read1, read2) is False):
            read1.set_tag("MP", "N")
            read2.set_tag("MP", "N")
            outHandle.write(read1)
            outHandle.write(read2)
            print("I am continuing over these records because the reads didn't overlap.")
            continue
        print("Making Layout1")
        Layout1 = Layout(read1)
        print("Making Layout2")
        Layout2 = Layout(read2)
        print("Merging Layouts")
        if(read1.pos < read2.pos):
            newReadString = Layout1.__csam__(Layout2)
            print("Merging didn't break everything!")
            print(newReadString)
            if(Layout1.MergeSuccess):
                outHandleMerge.write(newReadString)
            else:
                print("Bam record merge failed, write original reads to file.")
                read1.set_tag("MP", "F")
                read2.set_tag("MP", "F")
                outHandle.write(read1)
                outHandle.write(read2)
        else:
            newReadString = Layout2.__csam__(Layout1)
            print("Merging didn't break everything!")
            if(Layout2.MergeSuccess):
                print("New SAM read: %s" % newReadString)
                outHandleMerge.write(newReadString)
                print("And neither did writing!")
            else:
                read1.set_tag("MP", "F")
                read2.set_tag("MP", "F")
                outHandle.write(read1)
                outHandle.write(read2)
    outHandle.close()
    inHandle.close()
    return outBAM


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


def testLayout(cystr inBAM, cystr outBAM=None):
    cdef AlignedSegment_t read
    cdef Layout_t layout
    if(outBAM is None):
        outBAM = TrimExt(inBAM) + ".testLayout.bam"
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    outHandle = pysam.AlignmentFile(outBAM, "wb", template=inHandle)
    for read in inHandle:
        layout = Layout(read)
        read = layout.__read__()
        outHandle.write(read)
    inHandle.close()
    outHandle.close()
    return outBAM


cdef list FlattenCigar(list cigar):
    cdef size_t k
    cdef list retList = []
    for k, g in groupby(cigar, oig0):
        glist = list(g)
        retList.append((k, sum(map(oig1, glist))))
    return retList


cpdef cystr FlattenCigarString(cystr cigar):
    return cFlattenCigarString(cigar)


cdef cystr cFlattenCigarString(cystr cigar):
    """
    Flattens cigar strings for cases where it is needed.
    """
    cdef size_t count = 0
    cdef int runSum = 0
    cdef int tmpInt
    cdef char cigarOpChar, tmpCigarOpChar
    cdef cystr retStr, tmpIntStr
    cdef int cigarOpLen
    cdef int newCigarOpLen
    cdef py_array chars = cs_to_ia("".join(rsplit("[0-9]+", cigar)[1:]))
    cdef py_array Lengths = array('i', [int(tmpIntStr) for tmpIntStr in
                                        rsplit("[A-Z]", cigar)[:len(chars)]])
    cdef list outTupleList = []
    '''
    cdef py_array Lengths = array(
        "i", )
    '''
    for CigarOpChar, entries in groupby(zip(chars, Lengths), key=oig0):
        for tmpCigarOpChar, tmpInt in list(entries):
            runSum += tmpInt
        outTupleList.append(opLenToStr(CigarOpChar, runSum))
        runSum = 0
    return "".join(outTupleList)
