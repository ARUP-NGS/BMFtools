from MawCluster.SVUtils import *
from utilBMF import HTSUtils

import pysam
import numpy as np

import uuid


def XLocIntrachromosomalFusionCaller(inBAM,
                                     minMQ=0,
                                     minBQ=0,
                                     bedfile="default",
                                     minClustDepth=3,
                                     minPileupLen=20,
                                     outfile="default",
                                     ref="default"):
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
    if(bedfile == "default"):
        FacePalm("Capture bed file required for translocation detection.")
    if(ref == "default"):
        ThisIsMadness("Path to reference index required!")
    # Do calls for LI first.
    # Now looking for intrachromosomal translocations
    LIBamRecords = HTSUtils.LoadReadPairsFromFile(inBAM, SVTag="LI,ORB",
                                                  minMQ=minMQ, minBQ=minBQ)
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    header = inHandle.header
    print("Number of records meeting requirements: {}".format(
        len(LIBamRecords)))
    ContigList = list(set([pair.read1_contig for pair in LIBamRecords]))
    PutativeXLocs = []
    parsedBedfile = HTSUtils.ParseBed(bedfile)
    for contig in ContigList:
        print("Beginning contig: {}".format(contig))
        WorkingPairSet = [pair for pair in LIBamRecords
                          if pair.read1_contig == contig]
        Clusters = ClusterByInsertSize(WorkingPairSet)
        if(len(Clusters) == 0):
            continue
        PutXIntervals = PileupISClustersByPos(Clusters,
                                              minClustDepth=minClustDepth,
                                              bedfile=parsedBedfile,
                                              minPileupLen=minPileupLen,
                                              header=header)
        print("Number of putative events to check: {}".format(
            len(PutXIntervals)))
        print(repr(PutXIntervals))
        PutTransReadPairSets = [SVSupportingReadPairs(interval, inHandle)
                                for interval in PutXIntervals]
        for event in PutTransReadPairSets:
            [bedIntervalList,
             MeanDORList] = HTSUtils.CreateIntervalsFromCounter(
                HTSUtils.ReadPairListToCovCounter(event,
                                                  minClustDepth=minClustDepth,
                                                  minPileupLen=minPileupLen),
                minPileupLen=minPileupLen,
                contig=contig,
                bedIntervals=parsedBedfile)
            PutativeXLocs.append(PutativeXLoc(DORList=MeanDORList,
                                              intervalList=bedIntervalList,
                                              ReadPairs=event,
                                              bedIntervals=parsedBedfile,
                                              header=header,
                                              TransType="Intrachromosoma"
                                              "lRearrangement",
                                              inBAM=inBAM))
    XLocVCFLines = list(set([TranslocationVCFLine(xLoc, inBAM=xLoc.inBAM,
                                                  ref=ref)
                             for xLoc in PutativeXLocs if
                             xLoc.nsegments != 0]))
    for line in XLocVCFLines:
        if(line.NumPartners != 0 and line.TDIST >= 50000):
            print("TDIST: {}".format(line.TDIST))
            outHandle.write(line.ToString() + "\n")
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
