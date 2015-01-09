from MawCluster.SNVUtils import *
import pysam

"""
Programs which write VCFs.
Currently: SNVCrawler.

In development: SV
"""


def SNVCrawler(inBAM,
               bed="default",
               minMQ=0,
               minBQ=0,
               OutVCF="default",
               MaxPValue=1e-15,
               keepConsensus=False,
               reference="default",
               reference_is_path=False,
               commandStr="default",
               fileFormat="default",
               FILTERTags="default",
               INFOTags="default",
               FORMATTags="default"
               ):
    if(isinstance(bed, str) and bed != "default"):
        bed = HTSUtils.ParseBed(bed)
    if(OutVCF == "default"):
        OutVCF = inBAM[0:-4] + ".bmf.vcf"
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    outHandle = open(OutVCF, "w")
    outHandle.write(GetVCFHeader(fileFormat=fileFormat,
                                 FILTERTags=FILTERTags,
                                 commandStr=commandStr,
                                 reference=reference,
                                 reference_is_path=False,
                                 header=inHandle.header,
                                 INFOTags=INFOTags,
                                 FORMATTags=FORMATTags))

    if(bed != "default"):
        for line in bed:
            puIterator = inHandle.pileup(line[0], line[1],
                                         max_depth=30000,
                                         multiple_iterators=True)
            while True:
                try:
                    PileupColumn = puIterator.next()
                    PC = PCInfo(PileupColumn, minMQ=minMQ, minBQ=minBQ)
                    # print(PC.toString())
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
    return OutVCF
