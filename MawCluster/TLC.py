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
    MDCBamRecords = HTSUtils.LoadReadPairsFromFile(inBAM, SVTag="MDC,ORB")
    return None
