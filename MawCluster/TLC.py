from MawCluster.SVUtils import SVParamDict
from utilBMF.HTSUtils import ParseBed, printlog as pl, FacePalm, DefaultSamHeader
from utilBMF import HTSUtils

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
            assert isinstance(bedIntervals[0], str) and isinstance(
                bedIntervals[1], int)
        except AssertionError:
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
                                           ["InBed=" + self.IntervalInBed]]])


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
                 bedIntervals="default"):
        try:
            assert (isinstance(DORList, list) and
                    isinstance(intervalList, list) and
                    isinstance(ReadPairs, list))
        except AssertionError:
            pl("DORList: {}".format(repr(DORList)))
            pl("intervalList: {}".format(repr(intervalList)))
            pl("ReadPairs: {}".format(repr(ReadPairs)))
            FacePalm("One of the given arguments is "
                     "not a list. Abort mission!")
        assert isinstance(DORList[0], float)
        try:
            assert (isinstance(intervalList[0][0], str) and
                    isinstance(intervalList[0][1], int))
        except TypeError:
            FacePalm("Interval list is not a list"
                     " of lists of [str, int, int]!")
        except AssertionError:
            FacePalm("Interval list not in appropriate format.")
        try:
            assert isinstance(ReadPairs[0], HTSUtils.ReadPair)
        except AssertionError:
            FacePalm("ReadPairs argument is not a list of ReadPair objects!")
        self.segments = [XLocSegment(DOR=dor, bedIntervals=bedIntervals,
                                     interval=interval)
                         for dor, interval in zip(DORList, intervalList)]
        self.ID = str(uuid.uuid4().get_hex().upper()[0:12])

    def ToString(self):
        string = ("@PutativeTranslocationID={}\tContig\tStart "
                  "[0-based]\tStop\tMean DOR\n".format(self.ID))
        for segment in self.segments:
            string += segment + "\n"

    def WriteReads(self, outBAM="default"):
        if(outBAM == "default"):
            outBAM = self.ID + ".xloc.bam"
        for ReadPair in self.ReadPairs:
            HTSUtils.WritePairToFile(ReadPair, outfile=outBAM)


def ClusterByInsertSize(ReadPairs, distance="default", bedfile="default",
                        insDistance="default"):
    # Check that it's a list of ReadPairs
    assert isinstance(ReadPairs[0], HTSUtils.ReadPair)
    # Assert that these are all mapped to the same contig
    assert len(list(set([pair.read1_contig for pair in ReadPairs]))) == 1
    if(distance == "default"):
        distance = ReadPairs[0].read1.query_length * 2
        pl("No distance provided - default of 2 * "
           "read length set: {}".format(distance))
    if(insDistance == "default"):
        insDistance = ReadPairs[0].read1.query_length * 2
        pl("No insert size distance provided - default of 2 * "
           "read length set: {}".format(insDistance))
    ClusterList = []
    for pair in ReadPairs:
        # Checks to see if it's already in a cluster. If so, skip.
        if(len(ClusterList) != 0):
            for cluster in ClusterList:
                if(pair in cluster):
                    continue
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
                    if(HTSUtils.ReadPairsWithinDistance(p, l)):
                        NewCS.append(l)
            Cluster = list(set(NewCS))
            if(len(Cluster) == clusterSize):
                break
        del NewCS
        del InsertSizeClusterSet
        if(clusterSize >= 3):
            ClusterList.append(Cluster)
            print(repr(ClusterList))
    return ClusterList


def PileupISClustersByPos(ClusterList, minClustDepth=3,
                          bedfile="default", minPileupLen=10):
    if(isinstance(bedfile, str)):
        bedfile = HTSUtils.ParseBed(bedfile)
    assert isinstance(ClusterList[0][0], HTSUtils.ReadPair)
    from collections import Counter
    for cluster in ClusterList:
        print("Length of cluster: {}".format(len(cluster)))
    for cluster in ClusterList:
        PutativeEvents = []
        for cs in cluster:
            # Find locations at which each cluster piles up
            posList = []
            posListDuplex = []
            potTransPos = {}
            print(repr(cs) + " is the repr of the cluster I'm looking at")
            for pair in cs:
                R1Pos = pair.read1.get_reference_positions()
                R2Pos = pair.read2.get_reference_positions()
                posList += R1Pos
                posList += R2Pos
                if([pos for pos in R1Pos if pos in R2Pos]):
                    posListDuplex.append(pos)
            PosCounts = Counter(posList)
            PosDuplexCounts = Counter(posListDuplex)
            # decrement the counts for each position to account for
            # both reads in a pair mapping to the same location.
            for key in PosDuplexCounts.keys():
                PosCounts[key] += -1 * PosDuplexCounts[key]
            for key in PosCounts.keys():
                potTransPos[key] = []
            PosCounts = dict([i for i in PosCounts.items()
                              if i[1] >= minClustDepth])
            # How should minPileupLen be chosen?
            if(len(PosCounts) < minPileupLen):
                continue
            bedIntervalList, MeanDORList = HTSUtils.CreateIntervalsFromCounter(
                PosCounts, minPileupLen=minPileupLen,
                contig=ClusterList[0][0].read1_contig)
            PutativeEvent = PutativeXLoc(DORList=MeanDORList,
                                         intervalList=bedIntervalList,
                                         ReadPairs=cluster, bedLines=bedfile)
            PutativeEvents.append(PutativeEvent)
    return PutativeEvents
    # TODO: star the ones which are outside the bedfile.


def XLocIntrachromosomalFusionCaller(inBAM,
                                     minMQ=0,
                                     minBQ=0,
                                     bedfile="default",
                                     minClustDepth=3,
                                     minPileupLen=10,
                                     outfile="default"):
    """
    Makes calls for translocations using SV Tags placed during SVUtils
    BAM must have these tags in order to find structural variants.
    Name Sorting required.
    Minimum BQ is not recommended for translocation calls.
    """

    if(outfile == "default"):
        outfile = inBAM[0:-4] + ".putativeSV.txt"
        print("Outfile: {}".format(outfile))
    outHandle = open(outfile, "w")
    print("{}".format(isinstance(outHandle, file)))
    """
    if(isinstance(bedfile, str) is True and bedfile != "default"):
        bedLines = ParseBed(bedfile)
    """
    # Do calls for LI first.
    # Now looking for intrachromosomal translocations
    LIBamRecords = HTSUtils.LoadReadPairsFromFile(inBAM, SVTag="LI,ORB",
                                                  minMQ=minMQ, minBQ=minBQ)
    print("Number of records meeting requirements: {}".format(
        len(LIBamRecords)))
    ContigList = list(set([pair.read1_contig for pair in LIBamRecords]))
    for contig in ContigList:
        WorkingPairSet = [pair for pair in LIBamRecords
                          if pair.read1_contig == contig]
        Clusters = ClusterByInsertSize(WorkingPairSet)
        PutativeEvents = PileupISClustersByPos(Clusters,
                                               minClustDepth=minClustDepth,
                                               bedfile=bedfile,
                                               minPileupLen=minPileupLen)
        for event in PutativeEvents:
            outHandle.write(event.ToString())
            event.WriteReads()
    outHandle.close()
    """
    Draft calls complete for intrachromosomal rearrangements.
    Now debugging, then expansion to interchromosomal rearrangements.
    Step 4: Try to create the consensus sequence using soft-clipped reads
    Step 5 (?): Create a variant graph using glia and verify the translocation.
    """
    # MDCBamRecords = HTSUtils.LoadReadPairsFromFile(inBAM, SVTag="MDC,ORB")
    # For MDC, Repeat, except that before doing the (within distance) filter,
    # get sets of reads which align to the same set of different contigs.
    return None
