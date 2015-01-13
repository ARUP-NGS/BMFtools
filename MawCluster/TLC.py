from MawCluster.SVUtils import SVParamDict
from utilBMF.HTSUtils import ParseBed

import pysam
from utilBMF.HTSUtils import FacePalm
from utilBMF import HTSUtils


def XLocCaller(inBAM,
               minMQ=0,
               minBQ=0,
               bedfile="default"):
    """
    Makes calls for translocations using SV Tags placed during SVUtils
    BAM must have these tags in order to find structural variants.
    Name Sorting required.
    # Large Insert first, then MDC.
    """
    keyList = SVParamDict.keys()

    """
    if(isinstance(bedfile, str) is True and bedfile != "default"):
        bedLines = ParseBed(bedfile)
    """
    # Do calls for LI first.
    # Now looking for intrachromosomal translocations
    LIBamRecords = HTSUtils.LoadReadPairsFromFile(inBAM, SVTag="LI,ORB")

    #Step 1: Go through all read pairs and find read pairs which are within 2 read lengths of the first pair (Maybe make that distance a bit more lax when we get into the millions?)
    #Step 2: Go through Step 1's results and keep only the ones which overlap with at least one other member of the group.
    #Step 3: See how much depth we can get for all positions there. If we have, say, 3 families supporting the translocation, maybe make it a call
    #Step 4: Try to create the consensus sequence using soft-clipped reads
    #Step 5 (?): Create a variant graph using glia and verify the translocation.
    MDCBamRecords = HTSUtils.LoadReadPairsFromFile(inBAM, SVTag="MDC,ORB")
    # For MDC, Repeat, except that before doing the (within distance) filter, get sets of reads which align to the same set of different contigs.
    return None
