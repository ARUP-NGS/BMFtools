import subprocess
import pysam
import math

from MawCluster.PileupUtils import *

"""
This module contains a variety of tools for calling variants.
Currently, it primarily works with SNPs primarily with experimental
features present for structural variants
TODO: Filter based on variants supported by reads going both ways.
TODO: Make calls for SNPs, not just reporting frequencies.
Reverse Strand Fraction - RSF
Both Strands Support Variant - BS
Fraction of unmerged reads supporting variant - TF
Allele Fraction - AF
DP - Depth (Merged)
"""


class VCFLine:

    def __init__(self, AltAggregateObject):
        if(isinstance(AltAggregateObject, AltAggregateInfo) is False):
            raise HTSUtils.ThisIsMadness("VCFLine requires an "
                                         "AltAggregateInfo for initialization")
        self.CHROM = AltAggregateObject.contig
        self.POS = AltAggregateObject.pos
        self.CONS = AltAggregateObject.consensus
        self.ALT = AltAggregateObject.ALT
        self.QUAL = AltAggregateObject.SumBQScore
        if(AltAggregateObject.BothStrandSupport is True):
            self.QUAL *= 10
            # This is terribly arbitrary...
        self.minMQ = AltAggregateInfo.minMQ
        self.minBQ = AltAggregateInfo.minBQ
        self.InfoFields = {"DP": AltAggregateObject.MergedReads,
                           "AF": AltAggregateObject.MergedReads
                           / float(AltAggregateObject.DOC),
                           "TF": AltAggregateObject.TotalReads
                           / float(AltAggregateObject.DOCTotal),
                           "BS": AltAggregateObject.BothStrandSupport,
                           "RSF": AltAggregateObject.ReverseMergedReads
                           / AltAggregateObject.MergedReads,
                           "QA": AltAggregateObject.SumBQScore}


def SNVCrawler(inBAM,
               bed="default",
               minMQ=0,
               minBQ=0,
               VariantCallingTsv="default",
               MaxPValue=1e-15
               ):
    if(isinstance(bed, str)):
        bed = HTSUtils.ParseBed(bed)
        if(VariantCallingTsv == "default"):
            VariantCallingTsv = inBAM[0:-4] + ".vc.tsv"
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    outHandle = open(VariantCallingTsv, "w")
    outHandle.write("#Contig\tPosition (0-based)\tConsensus\t"
                    "ALT\tQUAL\tFILTER\tINFO\n")
    phredSumMin = -10 * math.log10(MaxPValue)
    for line in bed:
        puIterator = inHandle.pileup(line[0], line[1], line[2],
                                     max_depth=30000)
        while True:
            variantPasses = False
            try:
                PC = PCInfo(puIterator.next(), minMQ=minMQ, minBQ=minBQ)
            except ValueError:
                pl(("Pysam sometimes runs into errors during iteration which"
                    " are not handled with any elegance. Continuing!"))
                continue
            except StopIteration:
                pl("Finished iterations.")
                break
            for alt in PC.AltAlleleData:
                if alt.ALT == alt.consensus:
                    continue
                if(alt.SumBQScore >= phredSumMin):
                    variantPasses = True
                if(alt.BothStrandSupport is True):
                    variantPasses = True
    outHandle.close()
    inHandle.close()
    return None
