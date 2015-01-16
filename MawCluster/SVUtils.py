import pysam
from utilBMF.HTSUtils import printlog as pl
from utilBMF.HTSUtils import ThisIsMadness

"""
TODO: Make calls based on an analysis of SBI tags
SBI tag subsets - MDC/ORB, LI/ORB
"""
#  Made a dictionary for all structural variant candidate types
#  such that cycling through the list will be easier.
#  Extra field provided for each.
SVTestDict = {}


def LI_SV_Tag_Condition(read1, read2, tag, extraField=50000):
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
    from utilBMF.HTSUtils import ReadContainedInBed as ORIB
    return (sum([ORIB(read1, bedRef=bedRef), ORIB(read2, bedRef=bedRef)]) == 1)

SVTestDict['ORB'] = ORB_SV_Tag_Condition


def SBI_SV_Tag_Condition(read1, read2, tag, extraField="default"):
    """
    Gets reads where only one pair mapped inside the bed file
    and the insert size is either above a threshold (default: 50000),
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
    from utilBMF.HTSUtils import ReadContainedInBed as ORIB
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
SVParamDict['LI'] = 50000
SVParamDict['ORB'] = "default"
SVParamDict['SBI'] = ["default", 50000]


def GetSVRelevantRecordsPaired(inbam, SVBam="default",
                               bedfile="default",
                               supplementary="default",
                               maxInsert=50000,
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
