from utilBMF.HTSUtils import printlog as pl
from utilBMF.HTSUtils import ThisIsMadness
from utilBMF.HTSUtils import FacePalm, PysamToChrDict
from utilBMF import HTSUtils

import pysam

import numpy as np
import uuid


class XLocSegment:
    """
    Contains evidence regarding a potential translocation. This includes the
    depth of supporting reads, the regions involved, and a list of those reads.
    TODO: Eventually remove the assertions for speed.
    """
    def __init__(self, interval="default", DOR="default",
                 bedIntervals="default"):
        try:
            assert isinstance(interval[0], str) and isinstance(
                interval[1], int)
        except AssertionError:
            FacePalm("interval must be in ParseBed output format!")
        try:
            assert isinstance(bedIntervals[0][0], str) and isinstance(
                bedIntervals[0][1], int)
        except AssertionError:
            print(repr(bedIntervals))
            FacePalm("bedIntervals must be in ParseBed output format!")
        self.IntervalInBed = HTSUtils.IntervalOverlapsBed(
            interval, bedIntervals)
        self.interval = interval
        self.DOR = DOR
        self.bedIntervals = bedIntervals
        self.ID = str(uuid.uuid4().get_hex().upper()[0:8])

    def ToString(self):
        return "\t".join([str(i) for i in [["SegmentID=" + self.ID] +
                                           self.interval + [self.DOR] +
                                           ["InBed={}".format(
                                            self.IntervalInBed)]]])


class PutativeXLoc:
    """
    Contains a list of XLocSegment objects and has a ToString method which
    produces a line with all regions (most likely only two) involved in the
    rearrangement, along with a randomly generated uuid. This uuid is the
    basename for the bam file containing the supporting reads.
    TODO: Eventually remove the assertions for speed.
    """
    def __init__(self, DORList="default",
                 intervalList="default",
                 ReadPairs="default",
                 bedIntervals="default",
                 header="default",
                 TransType="UnspecifiedSV",
                 inBAM="default"):
        if(inBAM == "default"):
            ThisIsMadness("input BAM path required for SV VCF "
                          "Writing required.")
        self.TransType = TransType
        assert isinstance(header, dict)
        self.segments = [XLocSegment(DOR=dor, bedIntervals=bedIntervals,
                                     interval=interval)
                         for dor, interval in zip(DORList, intervalList)]
        self.ID = self.TransType + str(uuid.uuid4().get_hex().upper()[0:12])
        self.ReadPairs = ReadPairs
        self.intervals = intervalList
        self.bed = bedIntervals
        self.inBAM = inBAM
        self.nsegments = len(self.segments)

    def ToString(self):
        string = ("@PutativeTranslocationID={}\tContig\tStart "
                  "[0-based]\tStop\tMean DOR\n".format(self.ID))
        for segment in self.segments:
            if(isinstance(segment.ToString(), str)):
                string += segment.ToString() + "\n"
        return string

    def WriteReads(self, outBAM="default", header="default"):
        if(outBAM == "default"):
            outBAM = self.ID + ".xloc.bam"
        assert isinstance(header, dict)
        outHandle = pysam.AlignmentFile(outBAM, "wb", header=header)
        for ReadPair in self.ReadPairs:
            HTSUtils.WritePairToHandle(ReadPair, handle=outHandle)


def ClusterByInsertSize(ReadPairs, distance="default",
                        insDistance="default"):
    """
    Takes a list of ReadPair objects and return a list of lists of ReadPair
    objects. Each list of ReadPair objects has been clustered by insert size.
    The difference between this function and ClusterByInsertSize is that this
    only expands clusters with the reads outside of the bed file.
    Bug fix, 1/20/15: The expanding cluster feature expands by clustered
    insert sizes, not by overlaps.
    """
    # Check that it's a list of ReadPairs
    assert isinstance(ReadPairs[0], HTSUtils.ReadPair)
    # Assert that these are all mapped to the same contig
    assert len(list(set([pair.read1_contig for pair in ReadPairs]))) == 1
    if(insDistance == "default"):
        insDistance = ReadPairs[0].read1.query_length * 3
        pl("No insert size distance provided - default of "
           " 3 * read length set: {}".format(insDistance))
    ClusterList = []
    for pair in ReadPairs:
        # Checks to see if it's already in a cluster. If so, skip.
        # if(len(ClusterList) != 0):
        #    for cluster in ClusterList:
        #        if(pair in cluster):
        #            continue
        # Expand cluster iteratively based on current members.
        InsertSizeClusterSet = [l for l in ReadPairs if
                                HTSUtils.ReadPairsInsertSizeWithinDistance(
                                    pair, l, distance=insDistance)]
        NewCS = []
        Cluster = []
        while True:
            clusterSize = len(Cluster)
            for p in InsertSizeClusterSet:
                for l in ReadPairs:
                    if(HTSUtils.ReadPairsInsertSizeWithinDistance(
                            p, l, distance=insDistance)):
                        NewCS.append(l)
            Cluster = list(set(NewCS))
            if(len(Cluster) == clusterSize):
                break
        del NewCS
        del InsertSizeClusterSet
        if(clusterSize >= 3):
            ClusterList.append(Cluster)
            # print(repr(ClusterList))
    ClusterList = [i for i in ClusterList if len(i) != 0]
    return ClusterList


def SVSupportingReadPairs(bedInterval, bamHandle="default"):
    """
    Takes a bedInterval (chrom [str], start [int], stop [int],
    0-based, closed-end notation). Returns a list of ReadPair
    objects which support the putative translocation.
    This is typically then provided to ReadPairListToCovCounter
    """
    assert (isinstance(bedInterval[0], str) and
            isinstance(bedInterval[1], int))
    try:
        assert isinstance(bamHandle, pysam.calignmentfile.AlignmentFile)
    except AssertionError:
        pl("Working under the assumption that a string was "
           "provided instead of a BAM Handle.")
        bamHandle = pysam.AlignmentFile(bamHandle, "rb")
    ReadMateOutBed = []
    samFetch = bamHandle.fetch(bedInterval[0], bedInterval[1])
    while True:
        try:
            read = samFetch.next()
            if(read.reference_start > bedInterval[2] or
               PysamToChrDict[read.reference_id] != bedInterval[0]):
                break
            ReadMateOutBed.append(samFetch.next())
        except StopIteration:
            break
    ReadMateInBed = []
    for read in ReadMateOutBed:
        try:
            ReadMateInBed.append(bamHandle.mate(read))
        except ValueError:
            pl("Read mate not included, as it is unmapped.")
            pass
    ReadPairs = []
    for out, inside in zip(ReadMateOutBed, ReadMateInBed):
        if(out.is_read1):
            ReadPairs.append(HTSUtils.ReadPair(out, inside))
        else:
            ReadPairs.append(HTSUtils.ReadPair(inside, out))
    # Changed the list of read pairs to a list of sets of readpairs, in
    # case both reads are out of the bed region, so as to not artificially
    # inflate the number of supporting read families.
    return list(set(ReadPairs))


def PileupISClustersByPos(ClusterList, minClustDepth=2,
                          bedfile="default", minPileupLen=8,
                          header="default"):
    """
    Takes a list of lists of ReadPair objects which have been clustered by
    insert size, creates a list of intervals outside the bed capture region.
    These are then fed to SVSupportingReadPairs.
    TODO: Make it expand/cluster only in the direction which is outside
    the BED file.
    """
    assert len(ClusterList) != 0
    import copy
    if(isinstance(bedfile, str)):
        bedpath = copy.copy(bedfile)
        bedfile = HTSUtils.ParseBed(bedfile)
        pl("Bedfile parsed! Path: {}".format(bedpath))
        del bedpath
    try:
        assert isinstance(ClusterList[0][0], HTSUtils.ReadPair)
    except IndexError:
        print(repr(ClusterList))
        raise ThisIsMadness("Something is wrong!!!")
    from collections import Counter
    for cluster in ClusterList:
        print("Length of cluster: {}".format(len(cluster)))
    PotTransIntervals = []
    for cluster in ClusterList:
        PosCounts = HTSUtils.ReadPairListToCovCounter(
            cluster, minClustDepth=minClustDepth, minPileupLen=minPileupLen)
        # Make a list of coordinates for investigating
        [bedIntervalList, MeanDORList] = HTSUtils.CreateIntervalsFromCounter(
            PosCounts, minPileupLen=minPileupLen,
            contig=ClusterList[0][0].read1_contig,
            bedIntervals=bedfile)
        # Grab each region which lies outside of the bed file.
        RegionsToPull = []
        for bedLine, mean in zip(bedIntervalList, MeanDORList):
            if(HTSUtils.IntervalOverlapsBed(bedLine, bedIntervals=bedfile)
               is False):
                RegionsToPull.append(bedLine)
            else:
                FacePalm("Something's not working as hoped - regions not"
                         " in bed should have been filtered out already.")
        PotTransIntervals += RegionsToPull
    return PotTransIntervals


class TranslocationVCFLine:
    def __init__(self, PutativeXLocObj, ref="default", inBAM="default"):
        assert isinstance(PutativeXLocObj, PutativeXLoc)
        if(ref == "default"):
            raise ThisIsMadness("Reference must be provided for "
                                "creating a VCFLine.")
        if(isinstance(inBAM, str)):
            inBAM = pysam.AlignmentFile(inBAM, "rb")
        elif(isinstance(inBAM, pysam.calignmentfile.AlignmentFile) is False):
            raise ThisIsMadness("A source BAM file required for VCFLine.")
        segmentsInBed = [segment for segment in PutativeXLocObj.segments if
                         HTSUtils.IntervalOverlapsBed(segment.interval,
                                                      segment.bedIntervals)]
        segmentLengths = [segment.interval[2] - segment.interval[1] for segment
                          in PutativeXLocObj.segments]
        segmentLengthsInBed = [segment.interval[2] - segment.interval[1] for
                               segment in segmentsInBed]
        if(len(segmentsInBed) == 0):
            self.VCFRecSegment = PutativeXLocObj.segments[
                np.argmax(segmentLengths)]
        else:
            self.VCFRecSegment = segmentsInBed[np.argmax(segmentLengthsInBed)]
        self.partnerSegments = [segment for segment in
                                PutativeXLocObj.segments if
                                segment != self.VCFRecSegment]
        self.NumPartners = len(self.partnerSegments)
        self.CHROM = self.VCFRecSegment.interval[0]
        self.POS = self.VCFRecSegment.interval[1] + 1
        self.REF = "N"
        self.ALT = "<TRA>"
        self.QUAL = "QUALITY_IN_PROGRESS"
        self.ID = PutativeXLocObj.ID
        self.FILTER = "FILTER_IN_PROGRESS"
        # Number Merged Family Pairs supporting Structural Variant
        # Number Total Read Pairs Supporting Structural Variant
        StartsAndEnds = []
        for readpair in PutativeXLocObj.ReadPairs:
            StartsAndEnds += [readpair.read1.reference_start,
                              readpair.read2.reference_start,
                              readpair.read1.reference_end,
                              readpair.read1.reference_end]
        self.TDIST = sorted(StartsAndEnds)[-1] - sorted(StartsAndEnds)[0]
        self.InfoFields = {"NMFPSSV": len(PutativeXLocObj.ReadPairs),
                           "TYPE": "TRA",
                           "NTRPSSV": sum([int(pair.read1.opt("FM")) for pair
                                           in PutativeXLocObj.ReadPairs]),
                           "TDIST": self.TDIST}
        self.InfoFields["SVSEGS"] = ""
        for seg in self.partnerSegments:
            self.InfoFields["SVSEGS"] = "|".join(
                [str(i) for i in
                 [seg.interval[0], seg.interval[1],
                  seg.interval[2], seg.DOR]]) + ":"
        self.FormatFields = {}
        self.InfoStr = ";".join(
            ["=".join([key, str(self.InfoFields[key])])
             for key in sorted(self.InfoFields.keys())])
        if(len(self.FormatFields.keys()) == 0):
            self.FormatStr = "\t"
        else:
            self.FormatStr = (
                ":".join(sorted(self.FormatFields.keys())) +
                "\t" + ":".join(str(
                    self.FormatFields[key]) for key in sorted(
                        self.FormatFields.keys())))
        self.str = "\t".join([str(i) for i in [self.CHROM,
                                               self.POS,
                                               self.ID,
                                               self.REF,
                                               self.ALT,
                                               self.QUAL,
                                               self.FILTER,
                                               self.InfoStr,
                                               self.FormatStr]])

    def update(self):
        self.FormatKey = ":".join(sorted(self.FormatFields.keys()))
        self.FormatValue = ":".join([str(self.FormatFields[key])
                                     for key in
                                     sorted(self.FormatFields.keys())])
        self.FormatStr = (":".join(sorted(self.FormatFields.keys())) + "\t" +
                          ":".join(
                              str(self.FormatFields[key])
                              for key in sorted(self.FormatFields.keys())))
        self.InfoStr = ";".join([key + "=" + str(self.InfoFields[key])
                                 for key in sorted(self.InfoFields.keys())])

    def ToString(self):
        self.update()
        self.str = "\t".join([str(i) for i in [self.CHROM,
                                               self.POS,
                                               self.ID,
                                               self.REF,
                                               self.ALT,
                                               self.QUAL,
                                               self.FILTER,
                                               self.InfoStr,
                                               self.FormatStr]])
        return self.str

#  Made a dictionary for all structural variant candidate types
#  such that cycling through the list will be easier.
#  Extra field provided for each.
SVTestDict = {}


def LI_SV_Tag_Condition(read1, read2, tag, extraField=100000):
    maxInsert = int(extraField)
    return (abs(read1.tlen) >= maxInsert or abs(read2.tlen) >= maxInsert)

SVTestDict['LI'] = LI_SV_Tag_Condition


def MDC_SV_Tag_Condition(read1, read2, tag, extraField="default"):
    return (read1.reference_id != read2.reference_id)

SVTestDict['MDC'] = MDC_SV_Tag_Condition


def ORU_SV_Tag_Condition(read1, read2, tag, extraField="default"):
    return (sum([read1.is_unmapped, read2.is_unmapped]) == 1)

SVTestDict['ORU'] = ORU_SV_Tag_Condition


def MSS_SV_Tag_Condition(read1, read2, tag, extraField="default"):
    if(read1.reference_id == read2.reference_id):
        return ((sum([read1.is_reverse, read2.is_reverse]) != 1 and
                 read1.reference_id == read2.reference_id))
    else:
        return False

SVTestDict['MSS'] = MSS_SV_Tag_Condition


def ORB_SV_Tag_Condition(read1, read2, tag, extraField="default"):
    bedRef = extraField
    if(bedRef == "default"):
        raise ThisIsMadness("bedRef must be provded to run this test!")
    from utilBMF.HTSUtils import ReadOverlapsBed as ORIB
    return (sum([ORIB(read1, bedRef=bedRef), ORIB(read2, bedRef=bedRef)]) == 1)

SVTestDict['ORB'] = ORB_SV_Tag_Condition


def SBI_SV_Tag_Condition(read1, read2, tag, extraField="default"):
    """
    Gets reads where only one pair mapped inside the bed file
    and the insert size is either above a threshold (default: 100000),
    the reads are mapped to different contigs, or the reads are mapped
    to the same strand.
    extraField should contain a bedRef as formatted for ORB as field 0,
    and field 1 an integer for the minimum insert size to be marked
    as SBI.
    """
    try:
        SVTags = read1.opt("SV").split(',')
        if("ORB" in SVTags and ("LI" in SVTags or "MDC" in SVTags or
                                "MSS" in SVTags)):
            return True
    except KeyError:
        pass
    bedRef = extraField[0]
    if(bedRef == "default"):
        raise ThisIsMadness("bedRef must be provded to run this test!")
    from utilBMF.HTSUtils import ReadOverlapsBed as ORIB
    if(sum([ORIB(read1, bedRef=bedRef), ORIB(read2, bedRef=bedRef)]) == 1):
        if(read1.reference_id != read2.reference_id):
            return True
        if(abs(read1.tlen) > int(extraField[1])):
            return True
        if(sum([read1.is_reverse, read2.is_reverse]) != 1):
            return True
        return False
    return False

SVTestDict['SBI'] = SBI_SV_Tag_Condition
SVParamDict = {}
for key in SVTestDict.keys():
    SVParamDict[key] = ''
SVParamDict['LI'] = 100000
SVParamDict['ORB'] = "default"
SVParamDict['SBI'] = ["default", 100000]


def GetSVRelevantRecordsPaired(inbam, SVBam="default",
                               bedfile="default",
                               supplementary="default",
                               maxInsert=100000,
                               tempBAMPrefix="default",
                               FullBam="default",
                               summary="default"):
    """
    Requires a name-sorted, paired-end bam file where pairs have been kept
    together. (If a read is to be removed, its mate must also be removed.)
    Optionally, a supplementary bam file can be provided for additional
    information.
    Additionally, adds tags for different characteristics relevant
    to structural variants.
    If tempBAMPrefix is set, then reads relevant to each feature will be
    written to BAM files labeled accordingly.
    "SV" is the tag. It can hold multiple values, separated by commas.
    LI for Large Insert
    MDC for Mapped to Different Contig
    ORU for One Read Unmapped
    MSS for Mapped to Same Strand
    ORB for One Read In Bed Region
    SBI for having ORB and one of either MDC or LI
    (Spanning Bed with Improper pair)
    NF for None Found
    """
    if(SVBam == "default"):
        SVBam = '.'.join(inbam.split('.')[0:-1]) + '.sv.bam'
    if(FullBam == "default"):
        FullBam = '.'.join(inbam.split('.')[0:-1]) + '.SVmarked.bam'
    from utilBMF.HTSUtils import ParseBed
    bed = ParseBed(bedfile)
    global SVParamDict
    global SVTestDict
    SVParamDict['ORB'] = bed
    SVParamDict['LI'] = maxInsert
    SVParamDict['SBI'] = [bed, maxInsert]
    SVCountDict = {}
    for key in SVParamDict.keys():
        SVCountDict[key] = 0
    SVCountDict['NOSVR'] = 0  # "No Structural Variant Relevance"
    SVCountDict['SVR'] = 0  # "Structural Variant-Relevant"
    inHandle = pysam.AlignmentFile(inbam, "rb")
    SVOutHandle = pysam.AlignmentFile(SVBam, "wb", template=inHandle)
    FullOutHandle = pysam.AlignmentFile(FullBam, "wb", template=inHandle)
    FeatureList = sorted(SVTestDict.keys())
    pl("FeatureList: {}".format(FeatureList))
    for key in FeatureList:
        if(key not in SVParamDict.keys()):
            SVTestDict[key] = ""
    for read in inHandle:
        WritePair = False
        if(read.is_read1 is True):
            read1 = read
            continue
        if(read.is_read2 is True):
            read2 = read
        assert read1.query_name == read2.query_name
        for key in FeatureList:
            if(SVTestDict[key](
                    read1, read2, key, extraField=SVParamDict[key]) is True):
                try:
                    read1.setTag("SV", read1.opt("SV") + "," + key)
                    read2.setTag("SV", read2.opt("SV") + "," + key)
                    if("NF" in read1.opt("SV").split(",")):
                        read1.setTag(
                            "SV", ','.join([
                                i for i in read1.opt(
                                    "SV").split(
                                        ",") if i != "NF"]))
                    if("NF" in read2.opt("SV").split(",")):
                        read2.setTag(
                            "SV", ','.join([
                                i for i in read2.opt(
                                    "SV").split(
                                        ",") if i != "NF"]))
                except KeyError:
                    read1.setTag("SV", key)
                    read2.setTag("SV", key)
                WritePair = True
                SVCountDict[key] += 1
        if(WritePair is True):
            SVOutHandle.write(read1)
            SVOutHandle.write(read2)
            SVCountDict["SVR"] += 1
        else:
            SVCountDict["NOSVR"] += 1
            read1.setTag("SV", "NF")
            read2.setTag("SV", "NF")
        FullOutHandle.write(read1)
        FullOutHandle.write(read2)
    inHandle.close()
    SVOutHandle.close()
    FullOutHandle.close()
    SVCountDict["TOTAL"] = SVCountDict["SVR"] + SVCountDict["NOSVR"]
    for key in SVCountDict.keys():
        pl("Number of reads marked with key {}: {}".format(
            key, SVCountDict[key]))
    if(summary != "default"):
        writeSum = open(summary, "w")
        writeSum.write("#Category\tCount\tFraction\n")
        for key in SVCountDict.keys():
            writeSum.write(
                "{}\t{}\t{}\n".format(key,
                                      SVCountDict[key],
                                      SVCountDict[key] / float(
                                          SVCountDict['TOTAL'])))
        writeSum.close()
    return SVBam, FullBam
