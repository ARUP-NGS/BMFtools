from MawCluster.SVUtils import SVParamDict
from utilBMF.HTSUtils import ParseBed, printlog as pl
from utilBMF.HTSUtils import FacePalm
from utilBMF import HTSUtils

import numpy as np


def ClusterByInsertSize(ReadPairs, distance="default", bedfile="default",
                        insDistance="default"):
    # Check that it's a list of ReadPairs - or at least that the first one is.
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
    # ClusterList will be populated with a list of lists. Each list in
    # ClusterList is a list of reads whose insert sizes are within distance
    # of each other.
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
    return ClusterList


def ClusterISClustersByPos(ClusterList, minClustDepth=3,
                           bedfile="default"):
    assert isinstance(ClusterList[0][0], HTSUtils.ReadPair)
    from collections import Counter
    for cluster in ClusterList:
        PutativeClusters = []
        for cs in cluster:
            # Find locations at which each cluster piles up
            posList = []
            posListDuplex = []
            potTransPos = {}
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
            WarrantsFurtherInvestigation = False
            for key in PosCounts.keys():
                if(PosCounts[key] >= minClustDepth):
                    potTransPos[PosCounts[key]].append(key)
                    WarrantsFurtherInvestigation = True
            if(WarrantsFurtherInvestigation is False):
                continue
            PosCounts = potTransPos
            # TODO: PosCounts (ordered integer list) to list of continuous
            # bed intervals for list(PosCounts), then return bed interval
            # and the average depth of read in a tuple
            bedIntervalList = [["1", 5400, 50000], ["1", 400, 500]]
            # TODO: For each such interval, calculate the MeanDOR
            MeanDORList = [2.3333, 4311.4]
            PutativeClusters.append({"MeanDOR": MeanDORList,
                                     "BEDIntervals": bedIntervalList})

    return PutativeClusters
    # TODO: star the ones which are outside the bedfile.


def XLocCaller(inBAM,
               minMQ=0,
               minBQ=0,
               bedfile="default"):
    """
    Makes calls for translocations using SV Tags placed during SVUtils
    BAM must have these tags in order to find structural variants.
    Name Sorting required.
    Minimum BQ is not recommended for translocation calls.
    """
    keyList = SVParamDict.keys()

    """
    if(isinstance(bedfile, str) is True and bedfile != "default"):
        bedLines = ParseBed(bedfile)
    """
    # Do calls for LI first.
    # Now looking for intrachromosomal translocations
    LIBamRecords = HTSUtils.LoadReadPairsFromFile(inBAM, SVTag="LI,ORB",
                                                  minMQ=minMQ, minBQ=minBQ)
    ContigList = list(set([pair.read1_contig for pair in LIBamRecords]))
    for contig in ContigList:
        WorkingPairSet = [pair for pair in LIBamRecords
                          if pair.read1_contig == contig]
        Clusters = ClusterByInsertSize(WorkingPairSet)
    ISClusters = []
    WorkingISCluster = []
    """
    Step 1: Go through all read pairs and find read pairs which are within
        2 read lengths of the first pair (Maybe make that distance a bit more
        lax when we get into the millions?)
    Step 2: Go through Step 1's results and keep only the ones which overlap
        with at least one other member of the group.
    Step 3: See how much depth we can get for all positions there. If we have,
        say, 3 families supporting the translocation, maybe make it a call
    Step 4: Try to create the consensus sequence using soft-clipped reads
    Step 5 (?): Create a variant graph using glia and verify the translocation.
    """
    MDCBamRecords = HTSUtils.LoadReadPairsFromFile(inBAM, SVTag="MDC,ORB")
    # For MDC, Repeat, except that before doing the (within distance) filter, get sets of reads which align to the same set of different contigs.
    return None
