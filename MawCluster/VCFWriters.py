import logging

import pysam

from MawCluster.SNVUtils import *
from MawCluster.PileupUtils import pPileupColumn


"""
Programs which write VCFs.
Currently: SNVCrawler.

In development: SV

Note/TODO: Use two main filters for SNVs (after a min Q30):
1. Min FA (number of family members agreeing) [int=3]
2. Min Fraction Agreement Within Family [float=0.6667]

Other possible, but less critical options.
1. No LI/MDC/MSS/ORU/ORS reads.
2. Requiring at least (1?) duplex, even for larger panels?
3. Min PVFrac (Tosses out a read for up to a significant fraction
of its positions if certain parts of the read have a lower PV)


cfDNA limits:
1. No LI/MDC/MSS/ORU/ORS reads. (Probably very important)
2. Require N duplex read pairs [int = 2 (?)]

More thoughts?
"""


def SNVCrawler(inBAM,
               bed="default",
               minMQ=0,
               minBQ=0,
               OutVCF="default",
               MaxPValue=1e-30,
               keepConsensus=False,
               reference="default",
               reference_is_path=False,
               commandStr="default",
               fileFormat="default",
               FILTERTags="default",
               INFOTags="default",
               FORMATTags="default",
               writeHeader=True,
               minFracAgreed=0.0, minFA=2):
    pl("Command to reproduce function call: "
       "SNVCrawler({}, bed=\"{}\"".format(inBAM, bed) +
       ", minMQ={}, minBQ={}, OutVCF".format(minMQ, minBQ) +
       "=\"{}\", MaxPValue={}".format(OutVCF, MaxPValue) +
       ",keepConsensus={}, reference=".format(keepConsensus) +
       "\"{}\", reference_is_path={}".format(reference, reference_is_path) +
       "commandStr=\"{}\", fileFormat=\"{}\"".format(commandStr, fileFormat) +
       ", FILTERTags=\"{}\", INFOTags=\"{}\"".format(FILTERTags, INFOTags) +
       ", FORMATTags=\"{}\", writeHeader={}".format(FORMATTags, writeHeader) +
       ", minFracAgreed={}, minFA={})".format(minFracAgreed, minFA))
    if(bed != "default"):
        pl("Bed file used: {}".format(bed))
    if(isinstance(bed, str) and bed != "default"):
        bed = HTSUtils.ParseBed(bed)
    if(OutVCF == "default"):
        OutVCF = inBAM[0:-4] + ".bmf.vcf"
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    outHandle = open(OutVCF, "w+")
    if(writeHeader):
        try:
            outHandle.write(GetVCFHeader(fileFormat=fileFormat,
                                         FILTERTags=FILTERTags,
                                         commandStr=commandStr,
                                         reference=reference,
                                         reference_is_path=False,
                                         header=inHandle.header,
                                         INFOTags=INFOTags,
                                         FORMATTags=FORMATTags))
        except ValueError:
            pl("Looks like the RG header wasn't parseable by pysam - that's u"
               "sually an artefact of the clash between pysam and GATK's ways"
               "of working with RG fields.", level=logging.DEBUG)
            outHandle.write(GetVCFHeader(fileFormat=fileFormat,
                                         FILTERTags=FILTERTags,
                                         commandStr=commandStr,
                                         reference=reference,
                                         reference_is_path=False,
                                         INFOTags=INFOTags,
                                         FORMATTags=FORMATTags))
    if(bed != "default"):
        for line in bed:
            pl("Combing through bed region {}".format(line),
               level=logging.DEBUG)
            puIterator = inHandle.pileup(line[0], line[1],
                                         max_depth=30000,
                                         multiple_iterators=True)
            while True:
                try:
                    PileupColumn = pPileupColumn(puIterator.next())
                except StopIteration:
                    pl("Finished iterations.")
                    break
                '''
                except ValueError:
                    pl(("Pysam sometimes runs into errors during iteration w"
                        "hich aren't handled elegantly. Continuing!"),
                       level=logging.DEBUG)
                    continue
                '''
                PC = PCInfo(PileupColumn, minMQ=minMQ, minBQ=minBQ)
                #  pl("Position for pileup (0-based): {}".format(PC.pos),
                #     level=logging.DEBUG)
                if(line[2] <= PC.pos):
                    break
                VCFLineString = VCFPos(PC, MaxPValue=MaxPValue,
                                       keepConsensus=keepConsensus,
                                       reference=reference,
                                       minFracAgreed=minFracAgreed,
                                       minFA=minFA).ToString()
                if(len(VCFLineString) != 0):
                    outHandle.write(VCFLineString + "\n")
    else:
        puIterator = inHandle.pileup(max_depth=30000, multiple_iterators=True)
        while True:
            try:
                # Last command - 0 means iterator was where it crashed.
                PCpysam = pPileupColumn(puIterator.next())
                # Last command - 0 means the PCInfo call was where it crashed
                PC = PCInfo(PCpysam, minMQ=minMQ, minBQ=minBQ)
            except ValueError:
                raise ThisIsMadness("Trying to figure out what's going on in "
                                    "pysam's iterations.")
            except StopIteration:
                break
            # TODO: Check to see if it speeds up to not assign and only write.
            VCFLineString = VCFPos(PC, MaxPValue=MaxPValue,
                                   keepConsensus=keepConsensus,
                                   reference=reference,
                                   minFracAgreed=minFracAgreed,
                                   minFA=minFA).ToString()
            if(len(VCFLineString) != 0):
                outHandle.write(VCFLineString + "\n")
    return OutVCF


# Trying to "parallelize" this...
# I'll get around to it later.
"""


def SNVMinion(inBAM,
              bed="default",
              minMQ=0,
              minBQ=0,
              VCFLines="default",
              MaxPValue=1e-30,
              keepConsensus=False,
              reference="default",
              reference_is_path=False,
              commandStr="default",
              fileFormat="default",
              FILTERTags="default",
              INFOTags="default",
              FORMATTags="default"):
    if(isinstance(bed, str) and bed != "default"):
        bed = HTSUtils.ParseBed(bed)
    if(VCFLines == "default"):
        VCFLines = inBAM[0:-4] + ".bmf.vcf"
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    outHandle = open(VCFLines, "w+")
    if(bed != "default"):
        for line in bed:
            puIterator = inHandle.pileup(line[0], line[1],
                                         max_depth=30000,
                                         multiple_iterators=True)
            while True:
                try:
                    PileupColumn = puIterator.next()
                    PC = PCInfo(PileupColumn, minMQ=minMQ, minBQ=minBQ)
                    # print(PC.ToString())
                except ValueError:
                    pl(("Pysam sometimes runs into errors during iteration w"
                        "hich are not handled with any elegance. Continuing!"))
                    continue
                except StopIteration:
                    pl("Finished iterations.")
                    break
                if(line[2] <= PC.pos):
                    break
                VCFLineString = VCFPos(PC, MaxPValue=MaxPValue,
                                       keepConsensus=keepConsensus,
                                       reference=reference
                                       ).ToString()
                if(len(VCFLineString) != 0):
                    outHandle.write(VCFLineString + "\n")
    else:
        puIterator = inHandle.pileup(max_depth=30000)
        while True:
            try:
                PC = PCInfo(puIterator.next(), minMQ=minMQ, minBQ=minBQ)
            except ValueError:
                pl(("Pysam sometimes runs into errors during iteration which"
                    " are not handled with any elegance. Continuing!"))
                continue
            except StopIteration:
                break
            # TODO: Check to see if it speeds up to not assign and only write.
            VCFLineString = VCFPos(PC, MaxPValue=MaxPValue,
                                   keepConsensus=keepConsensus,
                                   reference=reference).ToString()
            if(len(VCFLineString) != 0):
                outHandle.write(VCFLineString + "\n")
    return VCFLines


def CallSNVCrawler():
    pass


def SNVMaster(inBAM,
              bed="default",
              minMQ=0,
              minBQ=0,
              VCFLines="default",
              MaxPValue=1e-30,
              keepConsensus=False,
              reference="default",
              reference_is_path=False,
              commandStr="default",
              fileFormat="default",
              FILTERTags="default",
              INFOTags="default",
              FORMATTags="default",
              ByContig=True,
              children=2):
    from subprocess import Popen
    if(isinstance(bed, str) and bed != "default"):
        bed = HTSUtils.ParseBed(bed)
    if(VCFLines == "default"):
        VCFLines = inBAM[0:-4] + ".bmf.vcf"
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    outHandle = open(VCFLines, "w")
    outHandle.write(GetVCFHeader(fileFormat=fileFormat,
                                 FILTERTags=FILTERTags,
                                 commandStr=commandStr,
                                 reference=reference,
                                 reference_is_path=False,
                                 header=inHandle.header,
                                 INFOTags=INFOTags,
                                 FORMATTags=FORMATTags))
    if(ByContig):
        contigList = list(set([line[0] for line in bed]))
        jobList = []
        for thread in range(int(children)):
            jobList.append(CallSNVCrawler(inBAM,
                           bed="default",
                           minMQ=0,
                           minBQ=0,
                           VCFLines="default",
                           MaxPValue=1e-30,
                           keepConsensus=False,
                           reference="default",
                           reference_is_path=False,
                           commandStr="default",
                           fileFormat="default",
                           FILTERTags="default",
                           INFOTags="default",
                           FORMATTags="default"))
    pass
"""
