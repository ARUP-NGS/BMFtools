#!/usr/bin/env python

from utilBMF.HTSUtils import printlog as pl
from utilBMF.HTSUtils import PysamToChrDict
from utilBMF.ErrorHandling import ThisIsMadness

import pysam
from cytoolz import frequencies as cyfreq


def MarkReadPairPositions(inBAM, outBAM="default"):

    """
    Input: Name-sorted BAM, without supplementary or secondary alignments.
    Output: Name-sorted BAM with RP tags which have the starts of the reads.
    """
    pl("Beginning MarkReadPairPositions")
    if(outBAM == "default"):
        outBAM = ".".join(inBAM.split('.')[0:-1]) + ".RP.bam"
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    outHandle = pysam.AlignmentFile(outBAM, "wb", template=inHandle)
    contigSets = []
    RPSetList = []
    while True:
        try:
            while True:
                read1 = inHandle.next()
                if(read1.is_secondary or read1.flag >> 11 == 1):
                    read1 = inHandle.next()
                else:
                    break
            while True:
                read2 = inHandle.next()
                if(read2.is_secondary or read2.flag >> 11 == 1):
                    read2 = inHandle.next()
                else:
                    break
            assert read2.qname == read1.qname
            coorString = ",".join(sorted([":".join([PysamToChrDict[
                read1.reference_id], str(read1.pos)]), ":".join([
                    PysamToChrDict[read2.reference_id], str(read2.pos)])]))
            contigSetStr = ",".join(sorted(
                [PysamToChrDict[read1.reference_id],
                 PysamToChrDict[read2.reference_id]]))
            RPSetList.append(coorString)
            read1.setTag("RP", coorString, "Z")
            read2.setTag("RP", coorString, "Z")
            read1.setTag("CS", contigSetStr, "Z")
            read2.setTag("CS", contigSetStr, "Z")
            outHandle.write(read1)
            outHandle.write(read2)
            contigSets.append(contigSetStr)
        except StopIteration:
            break
        except AssertionError:
            raise ThisIsMadness("Input BAM is not name-sorted!")
    RPSetCounts = cyfreq(RPSetList)
    RPSCHandle = open(inBAM[0:-3] + "rpsc.table", "w")
    for key in RPSetCounts.keys():
        RPSCHandle.write("\t".join([str(RPSetCounts[key]), str(key)]) + "\n")
    RPSCHandle.close()
    contigSets = list(set(contigSets))
    inHandle.close()
    outHandle.close()
    returnDict = {}
    returnDict["contigSets"] = contigSets
    returnDict["outBAM"] = outBAM
    return returnDict


def SortBamByRPTag(inBAM, outBAM="default", contigSets="default"):
    """
    Loads reads from each contig pair, sorts them by the RP tag, and
    then writes them to the output BAM in order of the RP string.
    BAM must be tagged with both RP and CS tags.
    """
    if(outBAM == "default"):
        outBAM = ".".join(inBAM.split('.')[0:-1]) + ".RPSort.bam"
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    outHandle = pysam.AlignmentFile(outBAM, "wb", template=inHandle)
    pl("Now sorting BAM by RP Tags")
    writeCount = 0
    for cSet in contigSets:
        inHandle = pysam.AlignmentFile(inBAM, "rb")
        readCount = 0
        pl("Working with contig set: {}".format(cSet.split(',')))
        workingSetR1 = []
        workingSetR2 = []
        while True:
            try:
                read = inHandle.next()
                if(read.opt("CS") != cSet):
                    continue
                if(read.is_read1):
                    workingSetR1.append(read)
                else:
                    workingSetR2.append(read)
                    readCount += 1
                continue
            except StopIteration:
                break
            except KeyError:
                raise ThisIsMadness("CS Key not present! "
                                    "Use MarkReadPairPositions")
            if(readCount % 1000 == 0):
                print("Number of read pairs read into memory"
                      " for contigSet {}: {}".format(cSet, readCount))
        workingSetR1 = sorted(workingSetR1, key=lambda x: x.opt("RP"))
        workingSetR2 = sorted(workingSetR2, key=lambda x: x.opt("RP"))
        for r1, r2 in zip(workingSetR1, workingSetR2):
            outHandle.write(r1)
            outHandle.write(r2)
            writeCount += 1
        print("{} read pairs written so far.".format(writeCount))
    inHandle.close()
    outHandle.close()
    return outBAM
