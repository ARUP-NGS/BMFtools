# cython: boundscheck=False, wraparound=False
## Standard library imports
from itertools import groupby, izip
from array import array
from operator import attrgetter as oag
import logging

## BMFTools imports
from .HTSUtils import TrimExt, printlog as pl, BamTag
from .ErrorHandling import ImproperArgumentError, ThisIsMadness as Tim

cimport cython

## DEFINES
oagtag = oag("tag")


cdef class LayoutPos:
    """
    Holds one layout position - either part of the reference,
    part of the read, or both.

    All of these fields are chars/ints.
    pos:
        -1 if "I", -1 * read.pos if "S"
    readPos:
        -1 if "D"
    mergeAgreed:
        1 = Not merged
        0 = Merge did not agree
        2 = Merge Agreed
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
                 bint merged=False, char mergeAgreed=1):
        self.pos = pos
        self.readPos = readPos
        self.operation = operation
        self.base = base
        self.quality = quality if(self.base != 78) else 0  # 78 == "N"
        self.agreement = agreement
        self.merged = merged
        self.mergeAgreed = mergeAgreed

    cpdef bint ismapped(self):
        return self.operation == 77  # 77 == "M"

    cdef bint getMergeAgreed(self):
        return self.mergeAgreed == 2

    cdef bint getMergeSet(self):
        return self.mergeAgreed != 1

    def __str__(self):
        return "%s|%s|%s|%s|%s|%s|%s|%s" % (
            self.pos, self.readPos, chr(self.base), chr(self.operation),
            self.quality, self.agreement, self.merged, self.mergeAgreed)


cdef class OldLayout(object):
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

    cdef py_array getSeqArr(self):
        """Returns a character array of the base calls
        if the base calls aren't "D" (68)
        """
        cdef char i
        cdef LayoutPos_t pos
        return array('B', [pos.base for pos in self.positions if
                           pos.operation != 68 and
                           pos.agreement > -1])

    cpdef cystr getSeq(self):
        return self.getSeqArr().tostring()

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

    cpdef int getAlignmentStart(self):
        cdef LayoutPos_t i
        for i in self.positions:
            if(i.operation == 77):  # operation is "M"
                return i.pos

    cpdef py_array getAgreement(self):
        """cpdef wrapper of cGetAgreement
        """
        return self.cGetAgreement()

    cdef py_array cGetAgreement(self):
        cdef LayoutPos_t pos
        # Ask if >= 0. My tests say it's ~1% faster to ask (> -1) than (>= 0).
        return array('i', [pos.agreement for
                           pos in self.positions if pos.agreement > -1])

    cdef py_array cGetQual(self):
        cdef LayoutPos_t pos
        # Ask if >= 0. My tests say it's ~1% faster to ask (> -1) than (>= 0).
        # pos.base == 66 for "B", which is a blank spot.
        # quality is set to less than 0 for an "N" cigar operation.
        return array('l', [pos.quality for
                           pos in self.positions if pos.operation != 66])

    cpdef py_array getQual(self):
        return self.cGetQual()

    cdef cystr cGetQualString(self):
        cdef int i
        return "".join([ph2chrInline(i) for i in self.cGetQual()])

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
                         lpos.merged],
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
            self.tagDict["MP"] = BamTag("MP", "Z", "T")

    cpdef update_tags(self):
        self.update_tags_()

    def update(self):
        cdef LayoutPos_t pos
        cdef int count
        self.update_tags()
        if(self.isMerged):
            raise NotImplementedError(
                "self.tlen should be set not to the length of seqArr, but to "
                "the number of cigar operation positions/aligned pairs/layou"
                "tpos's. Change it!")
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

    @cython.returns(py_array)
    def getOperations(self):
        """
        In [15]: %timeit [chr(p.operation) for p in l1.positions]
        100000 loops, best of 3: 9.35 us per loop

        In [16]: %timeit map(chr, [p.operation for p in l1.positions])
        100000 loops, best of 3: 8.69 us per loop
        """
        cdef LayoutPos_t pos
        # Skip this operation if the operation is "N"
        return array("B", [pos.operation for pos in
                           self.positions if
                           pos.operation != 78]).tostring()

    cdef cystr cGetCigarString(self):
        return "".join([str(len(list(g))) + k for
                        k, g in groupby(self.getOperations())])

    cpdef cystr getCigarString(self):
        return self.cGetCigarString()

    @cython.returns(list)
    def get_tags(self, object oagtag=oagtag):
        self.update_tags()
        return sorted(self.tagDict.itervalues(), key=oagtag)

    def getFlag(self):
        self.update()
        return self.flag

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
                      self.cGetCigarString(), self.rnext, self.pnext + 1,
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
        newRead.reference_id = ChrToRefIDInline(self.contig)
        newRead.query_sequence = self.getSeq()
        newRead.query_qualities = [93 if(i > 92) else i for
                                   i in self.cGetQual()]
        newRead.reference_start = self.getAlignmentStart()
        newRead.cigarstring = self.getCigarString()
        newRead.tlen = self.tlen
        newRead.mapq = self.mapq
        newRead.tags = [(BT.tag, BT.value) for BT in self.tagDict.itervalues()]
        newRead.next_reference_id = ChrToRefIDInline(self.rnext)
        newRead.pnext = self.pnext
        return newRead


    def __init__(self, pysam.calignmentfile.AlignedSegment rec,
                 list layoutPositions, int firstMapped):
        cdef tuple tag
        self.mapq = rec.mapq
        self.positions = layoutPositions
        self.firstMapped = firstMapped
        self.InitPos = rec.pos
        self.Name = rec.query_name
        self.contig = PysamToChrInline(rec.reference_id)
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
        self.rnext = PysamToChrInline(rec.mrnm)
        self.pnext = rec.mpos
        self.tlen = rec.tlen
        self.isMerged = (rec.has_tag("MP") and rec.opt("MP") == "T")
        self.is_reverse = rec.is_reverse


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
        if((not read1.is_proper_pair) or
           read1.reference_id != read2.reference_id or
           read1.is_unmapped or
           read2.is_unmapped or
           abs(read1.tlen) >= rLen2):
            ohw(read1)
            ohw(read2)
            continue
        try:
            assert(read1.query_name == read2.query_name)
        except AssertionError:
            raise ImproperArgumentError(
                "Bam is either not name sorted or you are missing a read from"
                " the pair around read # %s in the bam." % count)
        Layout1 = OldLayout.fromread(read1)
        Layout2 = OldLayout.fromread(read2)
        retLayout = MergeLayoutsToLayout(Layout1, Layout2)
        if(retLayout is None):
            read1.setTag("MP", "F")
            read2.setTag("MP", "F")
            ohw(read1)
            ohw(read2)
            continue
        ohw(retLayout.__read__())
    outHandle.close()
    inHandle.close()
    return outBAM


cpdef Layout_t MergeLayoutsToLayout(Layout_t L1, Layout_t L2):
    """
    Warning: This modifies L1 in-place under the assumption
    that the input arguments are to be ignored. Do not attempt
    to use L1 as the original layout.
    """
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
    """
    This should ideally only be called on reads which are already
    established as members in proper pairs.
    """
    cdef int pos1, pos2, readlen, end1, end2
    cdef LayoutPos_t tmpPos
    pos1 = L1.getRefPosForFirstPos()
    pos2 = L2.getRefPosForFirstPos()
    end1 = L1.cGetLastRefPos()
    end2 = L2.cGetLastRefPos()
    if(pos1 > end2 or pos2 > end1):
        return False
    else:
        return True


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
    offset = L2.cGetRefPosForFirstPos() - L1.cGetRefPosForFirstPos()
    return ListBool((L1[:offset] + [cMergePositions(pos1, pos2) for
                                    pos1, pos2 in
                                    izip(L1[offset:], L2)] +
                     L2[len(L1) - offset:]), True)


cdef LayoutPos_t cMergePositions(LayoutPos_t pos1, LayoutPos_t pos2):
    """Merges two positions. Order does matter - pos1 overrides pos2 when
    pos2 is soft-clipped.
    """
    if(pos1.base == pos2.base):
        if(pos2.operation == pos1.operation):
            return LayoutPos(pos1.pos, pos1.readPos, pos1.base, pos1.operation,
                             pos1.quality + pos2.quality,
                             pos1.agreement + pos2.agreement, merged=True,
                             mergeAgreed=2)
        elif(pos2.operation == 83):  # if pos2.operation is "S"
            return LayoutPos(pos1.pos, pos1.readPos, pos1.base,
                             pos1.operation,
                             pos1.quality + pos2.quality,
                             pos1.agreement + pos2.agreement,
                             merged=True, mergeAgreed=2)
        else:
            return LayoutPos(pos1.pos, pos1.readPos, 78, pos1.operation,
                             -137, -137,
                             merged=True, mergeAgreed=0)
    elif(pos1.operation == pos2.operation):
        if(pos1.quality > pos2.quality):
            return LayoutPos(
                pos1.pos, pos1.readPos, pos1.base, pos1.operation,
                pos1.quality - pos2.quality, pos1.agreement,
                merged=True, mergeAgreed=0)
        else:
            return LayoutPos(
                pos1.pos, pos1.readPos, pos2.base, pos1.operation,
                pos2.quality - pos1.quality, pos2.agreement,
                merged=True, mergeAgreed=0)
    else:
        return LayoutPos(pos1.pos, pos1.readPos, 78, pos1.operation, -137,
                         -137, merged=True, mergeAgreed=0)


@cython.returns(Layout_t)
def makeLayout(pysam.calignmentfile.AlignedSegment rec):
    return OldLayout(makeLayoutTuple)


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


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int getLayoutLen(AlignedSegment_t read):
    return len(read.aligned_pairs)


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
    CigarChar = CigarOpToCigarChar(cigarOp)
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
    return [LayoutPos(pos=x1, readPos=x0, operation=77,
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
