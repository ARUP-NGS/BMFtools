import pysam

from MawCluster.PileupUtils import PCInfo, AlleleAggregateInfo
from utilBMF import HTSUtils
from utilBMF.HTSUtils import printlog as pl

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

    def __init__(self,
                 AltAggregateObject,
                 MaxPValue=1e-15,
                 ID=".",
                 DOCMerged="default",
                 DOCTotal="default",
                 TotalAlleleFracStr="default",
                 MergedAlleleFracStr="default",
                 TotalAlleleCountStr="default",
                 MergedAlleleCountStr="default"):
        if(isinstance(AltAggregateObject, AlleleAggregateInfo) is False):
            raise HTSUtils.ThisIsMadness("VCFLine requires an AlleleAgg"
                                         "regateInfo for initialization")
        if(DOCMerged == "default"):
            raise HTSUtils.ThisIsMadness("DOC (Merged) required!")
        if(DOCTotal == "default"):
            raise HTSUtils.ThisIsMadness("DOC (Total) required!")
        self.CHROM = AltAggregateObject.contig
        self.POS = AltAggregateObject.pos
        self.CONS = AltAggregateObject.consensus
        self.ALT = AltAggregateObject.ALT
        self.QUAL = AltAggregateObject.SumBQScore
        if(AltAggregateObject.BothStrandSupport is True):
            self.QUAL *= 100
            # This is terribly arbitrary...
        self.ID = ID
        if(float(MaxPValue) > 10 ** float(self.QUAL / -10)):
            self.FILTER = "PASS"
        else:
            self.FILTER = "LowQual"
        if(self.ALT == AltAggregateObject.consensus):
            self.FILTER = "CONSENSUS"
        self.InfoFields = {"AF": AltAggregateObject.MergedReads
                           / float(AltAggregateObject.DOC),
                           "TF": AltAggregateObject.TotalReads
                           / float(AltAggregateObject.DOCTotal),
                           "BS": AltAggregateObject.BothStrandSupport,
                           "RSF": AltAggregateObject.ReverseMergedReads
                           / AltAggregateObject.MergedReads,
                           "MQM": AltAggregateObject.AveMQ,
                           "MQB": AltAggregateObject.AveBQ,
                           "MMQ": AltAggregateObject.minMQ,
                           "MBQ": AltAggregateObject.minBQ,
                           "QA": AltAggregateObject.SumBQScore,
                           "NUMALT": AltAggregateObject.NUMALT,
                           "BQF": AltAggregateObject.FailedBQReads,
                           "MQF": AltAggregateObject.FailedMQReads,
                           "TYPE": "snp",
                           "MVQ": MaxPValue}
        if(TotalAlleleCountStr != "default"):
            self.InfoFields["TACS"] = TotalAlleleCountStr
        if(TotalAlleleFracStr != "default"):
            self.InfoFields["TAFS"] = TotalAlleleFracStr
        if(MergedAlleleCountStr != "default"):
            self.InfoFields["MACS"] = MergedAlleleCountStr
        if(MergedAlleleFracStr != "default"):
            self.InfoFields["MAFS"] = MergedAlleleFracStr
        self.InfoStr = ";".join(["=".join(key, str(self.InfoFields[key]))
                                 for key in self.InfoFields.keys()])
        self.FormatFields = {"DP": DOCMerged,
                             "DPA": AltAggregateObject.MergedReads,
                             "DPT": DOCTotal,
                             "QA": AltAggregateObject.SumBQScore}
        self.FormatStr = (":".join(self.FormatFields.keys()) + "\t" +
                          ":".join(str(self.FormatFields[key])
                                   for key in self.FormatFields.keys()))
        self.str = ToString(self)

        def update():
            self.FormatStr = (":".join(self.FormatFields.keys()) + "\t" +
                              ":".join(str(self.FormatFields[key])
                                       for key in self.FormatFields.keys()))
            self.InfoStr = ";".join(["=".join(key, str(self.InfoFields[key]))
                                    for key in self.InfoFields.keys()])

        def ToString():
            self.update()
            self.str = "\t".join([str(i) for i in [self.CHROM,
                                                   self.POS,
                                                   self.ID,
                                                   self.CONS,
                                                   self.ALT,
                                                   self.QUAL,
                                                   self.FILTER,
                                                   self.InfoStr,
                                                   self.FormatStr]])
            return self.str


class VCFPos:

    def __init__(self, PCInfoObject, MaxPValue=1e-15, keepConsensus=True):
        if(isinstance(PCInfoObject, PCInfo) is False):
            raise HTSUtils.ThisIsMadness("VCFPos requires an "
                                         "PCInfo for initialization")
        self.pos = PCInfoObject.pos
        self.minMQ = PCInfoObject.minMQ
        self.consensus = PCInfoObject.consensus
        self.TotalAlleleFracStr = ",".join(
            [key +
             "->" +
             str(PCInfoObject.TotalFracFormatDict[key][0:6])
             for key in PCInfoObject.TotalFracFormatDict.keys()])
        self.MergedAlleleFracStr = ",".join(
            [key +
             "->" +
             str(PCInfoObject.MergedFracFormatDict[key][0:6])
             for key in PCInfoObject.MergedFracFormatDict.keys()])
        self.TotalAlleleCountStr = ",".join(
            [key +
             "->" +
             str(PCInfoObject.TotalCountFormatDict[key])
             for key in PCInfoObject.TotalCountFormatDict.keys()])
        self.MergedAlleleCountStr = ",".join(
            [key +
             "->" +
             str(PCInfoObject.MergedCountFormatDict[key])
             for key in PCInfoObject.MergedCountFormatDict.keys()])
        # Create VCFLines object using the calculated statistics
        self.VCFLines = [VCFLine(
            alt, TotalAlleleCountStr=self.TotalAlleleCountStr,
            MergedAlleleCountStr=self.MergedAlleleCountStr,
            TotalAlleleFracStr=self.TotalAlleleFracStr,
            MergedAlleleFracStr=self.MergedAlleleFracStr,
            DOCMerged=PCInfoObject.MergedReads,
            DOCTotal=PCInfoObject.TotalReads,
            MaxPValue=MaxPValue)
            for alt in PCInfoObject.AltAlleleData]
        self.str = self.ToString()
        self.keepConsensus = keepConsensus

        def ToString(self):
            [line.update() for line in self.VCFLines]
            if(self.keepConsensus is True):
                self.str = "\n".join([line.ToString()
                                      for line in
                                      self.VCFLines])
            else:
                self.str = "\n".join([line.ToString()
                                      for line in
                                      self.VCFLines
                                      if line.FILTER != "CONSENSUS"])
            return self.str


def SNVCrawler(inBAM,
               bed="default",
               minMQ=0,
               minBQ=0,
               VariantCallingTsv="default",
               MaxPValue=1e-15,
               keepConsensus=False
               ):
    if(isinstance(bed, str)):
        bed = HTSUtils.ParseBed(bed)
        if(VariantCallingTsv == "default"):
            VariantCallingTsv = inBAM[0:-4] + ".vc.tsv"
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    outHandle = open(VariantCallingTsv, "w")
    outHandle.write("##HEADER___NOT_WRITTEN")
    for line in bed:
        puIterator = inHandle.pileup(line[0], line[1], line[2],
                                     max_depth=30000)
        while True:
            try:
                PC = PCInfo(puIterator.next(), minMQ=minMQ, minBQ=minBQ)
            except ValueError:
                pl(("Pysam sometimes runs into errors during iteration which"
                    " are not handled with any elegance. Continuing!"))
                continue
            except StopIteration:
                pl("Finished iterations.")
                break
            VCFLineSet = VCFPos(PC, MaxPValue=MaxPValue,
                                keepConsensus=keepConsensus)
            # TODO: Check to see if it speeds up to not assign and only write.
            outHandle.write(VCFLineSet.ToString() + "\n")
    outHandle.close()
    inHandle.close()
    return None
