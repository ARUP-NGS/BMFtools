import subprocess
import decimal
import numpy as np

import pysam

from BMFUtils.HTSUtils import ThisIsMadness, printlog as pl
from BMFUtils import HTSUtils
import BMFUtils

'''
TODO: Filter based on variants supported by reads going both ways.
TODO: Create a tsv file containing all passing variants, frequencies.
TODO: Run a variant-caller on the sample as a whole to make a list
    from which to subtract variants and heterozygotes

'''


class VCFFile:

    """A simple VCFFile object, consisting of a header, a name for the file
    from which they came, and a list of all VCFRecords.
    Access header through self.header (list of strings, one for each line)
    Access VCF Entries through self.Records (list of VCFRecord objects)
    """

    def __init__(self, VCFEntries, VCFHeader, inputVCFName):
        self.sampleName = inputVCFName
        self.header = VCFHeader
        self.Records = VCFEntries
        self.sampleNamesArray = [inputVCFName]
        self.numSamples = len(self.sampleNamesArray)

    def cleanRecords(self):
        NewRecordList = [entry for entry in self.Records if entry.ALT != "X"]
        self.Records = NewRecordList

    def filter(self, filterOpt="default", param="default"):
        if(filterOpt == "default"):
            try:
                raise ValueError(
                    "Filter method required")
            except ValueError:
                print("Returning nothing.")
                return
        NewVCFEntries = []
        for entry in self.Records:
            if(VCFRecordTest(entry, filterOpt, param=param) is True):
                NewVCFEntries.append(entry)
        if(filterOpt == "bed"):
            if(param == "default"):
                try:
                    raise ValueError("Bed file required for bedfilter.")
                except ValueError:
                    print("Returning nothing.")
            filterOpt = filterOpt + "&" + param
        NewVCFFile = VCFFile(NewVCFEntries, self.header, self.sampleName
                             + "FilteredBy{}".format(filterOpt))
        # TODO: make a new VCFFile object based on location.
        return NewVCFFile

    def update(self):
        SetNames = []
        for record in self.Records:
            record.update()
            SetNames.append(record.VCFFilename)
        self.sampleNamesArray = list(set(SetNames))
        self.numSamples = len(self.sampleNamesArray)

    def len(self):
        try:
            return len(self.Records)
        except AttributeError:
            raise AttributeError("VCFFile object not initialized.")

    def write(self, Filename):
        FileHandle = open(Filename, "w")
        for headerLine in self.header:
            FileHandle.write("{}\n".format(headerLine))
        for VCFEntry in self.Records:
            FileHandle.write("{}\n".format(VCFEntry.toString()))
        FileHandle.close()
        return

    def GetLoFreqVariants(self, outVCF="default", replace=False):
        NewVCFEntries = [entry for entry in self.Records if VCFRecordTest(
                         entry, filterOpt="I16", param="3")]
        NewVCF = VCFFile(NewVCFEntries,
                         self.header,
                         self.sampleName.replace(".vcf", "") + ".lofreq.vcf")
        if(replace is True):
            self = NewVCF
        if(outVCF == "default"):
            pl("outVCF not set, not writing to file.")
        else:
            self.write(outVCF)
        return NewVCF


class VCFRecord:

    """A simple VCFRecord object, taken from an item from
    the list which ParseVCF returns as "VCFEntries" """

    def __init__(self, VCFEntry, VCFFilename):
        self.CHROM = VCFEntry[0]
        self.POS = VCFEntry[1]
        self.ID = VCFEntry[2]
        self.REF = VCFEntry[3]
        if("<X>" in VCFEntry[4]):
            self.ALT = VCFEntry[4].replace(",<X>", "")
        else:
            self.ALT = VCFEntry[4]
        '''
        if("<X>" != VCFEntry[4]):
            self.ALT = ','.join(VCFEntry[4].split(',').remove("<X>"))
        else:
            self.ALT = "<X>"
        '''
        self.QUAL = VCFEntry[5]
        self.FILTER = VCFEntry[6]
        self.INFO = VCFEntry[7]
        self.InfoKeys = [entry.split(
            '=')[0] for entry in self.INFO.split(';') if entry != "NoCG"]
        self.InfoValues = []
        self.InfoUnpaired = []
        for entry in self.INFO.split(';'):
            #  print("entry: {}. INFO: {}".format(entry, self.INFO))
            try:
                self.InfoValues.append(entry.split('=')[1])
            except IndexError:
                self.InfoUnpaired
                continue
        #  print(self.InfoValues)
        #  Might not reproduce the original information when written to file.
        tempValArrays = [entry.split(',') for entry in self.InfoValues]
        try:
            self.InfoValArrays = [
                [entry for entry in array] for array in tempValArrays]
        except ValueError:
            self.InfoValArrays = [
                [entry for entry in array] for array in tempValArrays]
        self.InfoDict = dict(zip(self.InfoKeys, self.InfoValues))
        self.InfoArrayDict = dict(zip(self.InfoKeys, self.InfoValArrays))
        try:
            self.InfoArrayDict['I16'] = [
                int(i) for i in self.InfoArrayDict['I16']]
        except ValueError:
            if("I16" in self.InfoArrayDict.keys()):
                self.InfoArrayDict['I16'] = [
                    int(decimal.Decimal(i)) for i in self.InfoArrayDict['I16']]
        try:
            self.FORMAT = VCFEntry[8]
        except IndexError:
            self.FORMAT = ""
        try:
            self.GENOTYPE = VCFEntry[9]
        except IndexError:
            self.GENOTYPE = ""
        self.GenotypeDict = dict(
            zip(self.FORMAT.split(':'), self.GENOTYPE.split(':')))
        self.GenotypeKeys = self.FORMAT.split(':')
        self.GenotypeValues = self.GENOTYPE.split(':')
        self.Samples = [""]
        if(len(VCFEntry) > 10):
            for field in VCFEntry[10:]:
                self.Samples.append(field)
        self.VCFFilename = VCFFilename
        if(len(self.Samples) == 0):
            recordStr = '\t'.join([self.CHROM,
                                   self.POS,
                                   self.ID,
                                   self.REF,
                                   self.ALT,
                                   self.QUAL,
                                   self.FILTER,
                                   self.INFO,
                                   self.FORMAT,
                                   self.GENOTYPE])
        else:
            sampleStr = "\t".join(self.Samples)
            recordStr = '\t'.join([self.CHROM,
                                   self.POS,
                                   self.ID,
                                   self.REF,
                                   self.ALT,
                                   self.QUAL,
                                   self.FILTER,
                                   self.INFO,
                                   self.FORMAT,
                                   self.GENOTYPE,
                                   sampleStr])
        self.str = recordStr.strip()

    def update(self):
        self.InfoValues = [','.join(
            InfoValArray) for InfoValArray in self.InfoValArrays]
        infoEntryArray = [InfoKey + "=" + InfoValue for InfoKey,
                          InfoValue in zip(self.InfoKeys, self.InfoValues)]
        self.INFO = ';'.join(infoEntryArray) + ';'.join(self.InfoUnpaired)
        self.InfoKeys = [entry.split('=')[0] for entry in self.INFO.split(';')]
        self.InfoValues = [
            entry.split('=')[1] for entry in self.INFO.split(';')]
        tempValArrays = [entry.split(',') for entry in self.InfoValues]
        try:
            self.InfoValArrays = [
                [entry for entry in array] for array in tempValArrays]
        except ValueError:
            self.InfoValArrays = [
                [entry for entry in array] for array in tempValArrays]
        self.InfoDict = dict(zip(self.InfoKeys, self.InfoValues))
        self.InfoArrayDict = dict(zip(self.InfoKeys, self.InfoValArrays))
        if('I16' in self.InfoArrayDict.keys()):
            try:
                self.InfoArrayDict['I16'] = [
                    int(i) for i in self.InfoArrayDict['I16']]
            except ValueError:
                if("I16" in self.InfoArrayDict.keys()):
                    self.InfoArrayDict['I16'] = [
                        int(decimal.Decimal(
                            i)) for i in self.InfoArrayDict['I16']]
        self.GenotypeKeys = self.FORMAT.split(':')
        self.GenotypeValues = self.GENOTYPE.split(':')
        self.FORMAT = ":".join(self.GenotypeKeys)
        self.GENOTYPE = ":".join(self.GenotypeValues)
        self.GenotypeDict = dict(
            zip(self.FORMAT.split(':'), self.GENOTYPE.split(':')))
        if(len(self.Samples) == 0):
            recordStr = '\t'.join([self.CHROM, self.POS,
                                   self.ID, self.REF, self.ALT, self.QUAL,
                                   self.FILTER, self.INFO, self.FORMAT,
                                   self.GENOTYPE])
        else:
            sampleStr = "\t".join(self.Samples)
            recordStr = '\t'.join([self.CHROM, self.POS, self.ID,
                                   self.REF, self.ALT, self.QUAL,
                                   self.FILTER, self.INFO, self.FORMAT,
                                   self.GENOTYPE, sampleStr])
        self.str = recordStr.strip()

    def toString(self):
        self.update()
        return self.str


class PRInfo:
    '''
    Created from a pysam.PileupRead object.
    Holds family size, SV tags, base quality,
    mapping quality, and base.
    '''
    def __init__(self, PileupRead):
        try:
            self.FM = int(PileupRead.alignment.opt("FM"))
        except KeyError:
            self.FM = 1
        try:
            self.SVTags = PileupRead.alignment.opt("SV").split(",")
        except KeyError:
            self.SVTags = "NF"
            # print("Warning: SV Tags unset.")
        self.BaseCall = PileupRead.alignment.query_sequence[
            PileupRead.query_position]
        self.BQ = PileupRead.alignment.query_qualities[
            PileupRead.query_position]
        self.MQ = PileupRead.alignment.mapq
        self.is_reverse = PileupRead.alignment.is_reverse


def is_reverse_to_str(boolean):
    if(boolean is True):
        return "reverse"
    elif(boolean is False):
        return "forward"
    else:
        return "unmapped"


class AltAggregateInfo:
    '''
    Class which holds summary information for a given alt allele
    at a specific position. Meant to be used as part of the PCInfo
    class as values in the VariantDict.
    All alt alleles in this set of reads should be identical.
    recList must be a list of PRInfo objects.
    '''
    def __init__(self, recList, acceptMQ0=False, consensus="default",
                 mergedSize="default",
                 totalSize="default"):
        import collections
        if(consensus == "default"):
            raise ThisIsMadness("A consensus nucleotide must be provided.")
        if(mergedSize == "default"):
            raise ThisIsMadness(("mergedSize must be provided: number of "
                                 "PRInfo at given position."))
        if(totalSize == "default"):
            raise ThisIsMadness(("mergedSize must be provided: number of "
                                 "PRInfo at given position."))
        # Check that all alt alleles are identical
        self.recList = recList
        self.acceptMQ0 = acceptMQ0
        try:
            assert(sum([rec.BaseCall == recList[
                0].BaseCall for rec in recList]) == len(recList))
        except AssertionError:
            # print("recList repr: {}".format(repr(recList)))
            print("Alt alleles: {}".format([i.BaseCall for i in recList]))
            raise ThisIsMadness(
                "AltAggregateInfo requires that all alt alleles agree.")
        self.TotalReads = np.sum([rec.FM for rec in recList])
        self.MergedReads = len(recList)
        self.AveFamSize = float(self.TotalReads) / self.MergedReads
        if(acceptMQ0 is False):
            self.SumBQScore = sum([rec.BQ for rec in recList if rec.MQ != 0])
        else:
            self.SumBQScore = sum([rec.BQ for rec in recList])
        if(acceptMQ0 is False):
            self.SumMQScore = sum([rec.MQ for rec in recList if rec.MQ != 0])
        else:
            self.SumMQScore = sum([rec.MQ for rec in recList])
        self.AveMQ = float(self.SumMQScore) / len(self.recList)
        self.AveBQ = float(self.SumBQScore) / len(self.recList)
        self.ALT = recList[0].BaseCall
        self.consensus = consensus
        self.transition = "->".join([consensus, self.ALT])

        self.strandedTransitions = {}
        self.strandedTransitionDict = collections.Counter(
            ["&&".join([self.transition, is_reverse_to_str(
                rec.is_reverse)]) for rec in self.recList])
        self.strandedTotalTransitionDict = collections.Counter(
            ["&&".join([self.transition, is_reverse_to_str(
                rec.is_reverse)]) for rec in self.recList]*rec.FM)
        self.strandedMergedTransitionDict = collections.Counter(
            ["&&".join([self.transition, is_reverse_to_str(
                rec.is_reverse)]) for rec in self.recList])
        self.StrandCountsDict = {}
        self.StrandCountsDict['reverse'] = sum([
            rec.is_reverse for rec in self.recList])
        self.StrandCountsDict['forward'] = sum([
            rec.is_reverse is False for rec in self.recList])
        self.StrandCountsTotalDict = {}
        self.StrandCountsTotalDict['reverse'] = sum(
            [rec.FM for rec in self.recList if rec.is_reverse is True])
        self.StrandCountsTotalDict['forward'] = sum(
            [rec.FM for rec in self.recList if rec.is_reverse is False])
        self.TotalAlleleFrequency = self.TotalReads / totalSize
        self.MergedAlleleFrequency = self.MergedReads / mergedSize


class PCInfo:
    '''
    Takes a pysam.PileupColumn covering one base in the reference
    and makes a new class which has "reference" (inferred from
    consensus) and a list of PRData (one for each read).
    '''
    def __init__(self, PileupColumn, acceptMQ0=False):
        from collections import Counter
        PysamToChrDict = BMFUtils.HTSUtils.GetRefIdDicts()['idtochr']
        self.contig = PysamToChrDict[PileupColumn.reference_id]
        self.pos = PileupColumn.reference_pos
        self.Records = [PRInfo(
            pileupRead) for pileupRead in PileupColumn.pileups]
        self.MergedReads = len(self.Records)
        try:
            self.TotalReads = sum([rec.FM for rec in self.Records])
        except KeyError:
            self.TotalReads = self.MergedReads
        self.consensus = Counter(
            [rec.BaseCall for rec in self.Records]).most_common(1)[0][0]
        self.VariantDict = {}
        for alt in list(set([rec.BaseCall for rec in self.Records])):
            self.VariantDict[alt] = [
                rec for rec in self.Records if rec.BaseCall == alt]
        # for key in self.VariantDict.keys():
        #     print("Key: {}. Value: {}.".format(key, self.VariantDict[key]))
        self.AltAlleleData = [AltAggregateInfo(
                              self.VariantDict[key],
                              acceptMQ0=acceptMQ0,
                              consensus=self.consensus,
                              mergedSize=self.MergedReads,
                              totalSize=self.TotalReads
                              ) for key in self.VariantDict.keys()]
        self.TotalFracDict = {}
        for alt in self.AltAlleleData:
            self.TotalFracDict[
                alt.ALT] = float(alt.TotalReads) / self.TotalReads
        self.MergedFracDict = {}
        for alt in self.AltAlleleData:
            self.MergedFracDict[
                alt.ALT] = float(alt.MergedReads) / self.MergedReads
        TransMergedCounts = {}
        for alt in self.AltAlleleData:
            try:
                TransMergedCounts[alt.transition] += alt.MergedReads
            except KeyError:
                TransMergedCounts[alt.transition] = alt.MergedReads
        self.TransMergedCounts = TransMergedCounts
        TransTotalCounts = {}
        for alt in self.AltAlleleData:
            try:
                TransTotalCounts[alt.transition] += alt.TotalReads
            except KeyError:
                TransTotalCounts[alt.transition] = alt.TotalReads
        self.TransTotalCounts = TransTotalCounts
        self.StrandedTransTotalCounts = {}
        for alt in self.AltAlleleData:
                for trans in alt.strandedTotalTransitionDict.keys():
                    try:
                        self.StrandedTransTotalCounts[
                            trans] += alt.strandedTotalTransitionDict[trans]
                    except KeyError:
                        self.StrandedTransTotalCounts[
                            trans] = alt.strandedTotalTransitionDict[trans]
        self.StrandedTransMergedCounts = {}
        for alt in self.AltAlleleData:
                for trans in alt.strandedMergedTransitionDict.keys():
                    try:
                        self.StrandedTransMergedCounts[
                            trans] += alt.strandedMergedTransitionDict[trans]
                    except KeyError:
                        self.StrandedTransMergedCounts[
                            trans] = alt.strandedMergedTransitionDict[trans]
        self.MergedAlleleDict = {"A" : 0, "C" : 0, "G": 0, "T": 0}
        self.TotalAlleleDict = {"A" : 0, "C" : 0, "G": 0, "T": 0}
        for alt in self.AltAlleleData:
            self.MergedAlleleDict[alt.consensus] = alt.MergedReads
            self.TotalAlleleDict[alt.consensus] = alt.TotalReads
        self.MergedAlleleFreqDict = {"A" : 0., "C" : 0., "G": 0., "T": 0.}
        self.TotalAlleleFreqDict = {"A" : 0., "C" : 0., "G": 0., "T": 0.}
        for key in self.MergedAlleleDict:
            self.MergedAlleleFreqDict[key] = self.MergedAlleleDict[
                key] / float(self.MergedReads)
            self.TotalAlleleFreqDict[key] = self.TotalAlleleDict[
                key] / float(self.TotalReads)
        
                
    def toString(self, header=False):
        outStr = ""
        if(header is True):
            outStr = ("#Chr\tPos (0-based)\tRef (Consensus)\tAlt\tTotal "
                      "Reads\tMerged Reads\tTotal Allele Frequency\tMerged "
                      "Allele Frequency\tReverse Total Reads\tForward Total"
                      " Reads\tReverse Merged Reads\tForward Merged Reads"
                      "\tFraction Of Total Reads\t"
                      "Fraction Of Merged Reads\tAverage "
                      "Family Size\t"
                      "BQ Sum\tBQ Mean\tMQ Sum\tMQ Mean\n")
        for alt in self.AltAlleleData:
            outStr += '\t'.join([str(i) for i in [self.contig, self.pos,
                                 self.consensus,
                                 alt.ALT,
                                 alt.TotalReads,
                                 alt.MergedReads,
                                 alt.TotalAlleleFrequency,
                                 alt.MergedAlleleFrequency,
                                 alt.StrandCountsTotalDict['reverse'],
                                 alt.StrandCountsTotalDict['forward'],
                                 alt.StrandCountsDict['reverse'],
                                 alt.StrandCountsDict['forward'],
                                 self.TotalFracDict[alt.ALT],
                                 self.MergedFracDict[alt.ALT],
                                 alt.AveFamSize,
                                 alt.SumBQScore,
                                 alt.AveBQ,
                                 alt.SumMQScore,
                                 alt.AveMQ]]) + "\n"
        self.str = outStr
        return self.str


def CustomPileupFullGenome(inputBAM,
                           PileupTsv="default",
                           TransitionTable="default",
                           StrandedTTable="default",
                           progRepInterval=1000):
    '''
    A pileup tool for creating a tsv for each position in the bed file.
    Used for calling SNPs with high confidence.
    Also creates several tables:
    1. Counts for all consensus-->alt transitions (By Total and Merged reads)
    2. Counts for the above, specifying strandedness
    3. Number of Merged Reads supporting each allele
    '''
    TransTotalDict = {}
    TransMergedDict = {}
    StrandedTransTotalDict = {}
    StrandedTransMergedDict = {}
    NumTransitionsTotal = 0
    NumTransitionsMerged = 0
    MergedReadsProcessed = 0
    TotalReadsProcessed = 0
    import os.path
    if(PileupTsv == "default"):
        PileupTsv = inputBAM[0:-4] + '.Pileup.tsv'
    if(TransitionTable == "default"):
        TransitionTable = inputBAM[0:-4] + '.Trans.tsv'
    if(StrandedTTable == "default"):
        StrandedTTable = inputBAM[0:-4] + '.StrandedTrans.tsv'
    if(os.path.isfile(inputBAM + ".bai") is False):
        pl("No bam index found for input bam - creating!")
        try:
            subprocess.check_call(['samtools', 'index', inputBAM])
        except subprocess.CalledProcessError:
            pl("Couldn't index BAM - coor sorting, then indexing!")
            inputBAM = HTSUtils.CoorSortAndIndexBam(inputBAM, uuid=True)
    NumProcessed = 0  # Number of processed positions in pileup
    bamHandle = pysam.AlignmentFile(inputBAM, "rb")
    PileupHandle = open(PileupTsv, "w")
    TransHandle = open(TransitionTable, "w")
    StrandedTransHandle = open(StrandedTTable, "w")
    FirstLine = True
    for pileupColumn in bamHandle.pileup():
        NumProcessed += 1
        if((NumProcessed) % progRepInterval == 0):
            pl("Number of positions processed: {}".format(
                NumProcessed))
            pl("Total reads processed: {}".format(TotalReadsProcessed))
            pl("Merged reads processed: {}".format(MergedReadsProcessed))
        PColSum = PCInfo(pileupColumn)
        MergedReadsProcessed += PColSum.MergedReads
        TotalReadsProcessed += PColSum.TotalReads
        if(FirstLine is True):
            PileupHandle.write(PColSum.toString(header=True))
            FirstLine = False
        else:
            PileupHandle.write(PColSum.toString())
        for key in PColSum.TransMergedCounts.keys():
            try:
                TransMergedDict[
                    key] += PColSum.TransMergedCounts[key]
            except KeyError:
                TransMergedDict[
                    key] = PColSum.TransMergedCounts[key]
            NumTransitionsMerged += PColSum.TransMergedCounts[
                key]
        for key in PColSum.TransTotalCounts.keys():
            try:
                TransTotalDict[
                    key] += PColSum.TransTotalCounts[key]
            except KeyError:
                TransTotalDict[
                    key] = PColSum.TransTotalCounts[key]
            NumTransitionsTotal += PColSum.TransTotalCounts[
                key]
        for key in PColSum.StrandedTransMergedCounts.keys():
            try:
                StrandedTransMergedDict[
                    key] += PColSum.StrandedTransMergedCounts[
                        key]
            except KeyError:
                StrandedTransMergedDict[
                    key] = PColSum.StrandedTransMergedCounts[
                        key]
        for key in PColSum.StrandedTransTotalCounts.keys():
            try:
                StrandedTransTotalDict[
                    key] += PColSum.StrandedTransTotalCounts[
                        key]
            except KeyError:
                StrandedTransTotalDict[
                    key] = PColSum.StrandedTransTotalCounts[
                        key]
    TransHandle.write(("Transition\tTotal Reads With Transition (Unflattened)"
                       "\tMerged Reads With Transition\tFraction Of Total "
                       "Transitions\tFraction Of Merged Transitions\n"))
    for key in TransTotalDict.keys():
        TransHandle.write("{}\t{}\t{}\n".format(key,
                                                TransTotalDict[key],
                                                TransMergedDict[key],
                                                TransTotalDict[key]/float(
                                                    NumTransitionsTotal),
                                                TransMergedDict[key]/float(
                                                    NumTransitionsMerged)))
    StrandedTransHandle.write(("Transition+Strandedness\tTotal Reads "
                               "(Unflattened)\tMergedReads With Transition\t"
                               "Fraction Of Total (Unflattened) Transitions"
                               "\tFraction of Merged Transitions\n"))
    for key in StrandedTransTotalDict.keys():
        StrandedTransHandle.write("{}\t{}\t{}\n".format(key,
                                  StrandedTransTotalDict[key],
                                  StrandedTransMergedDict[key],
                                  StrandedTransTotalDict[key] / float(
                                      NumTransitionsTotal),
                                  StrandedTransMergedDict[key] / float(
                                      NumTransitionsMerged),
                                  ))
    pl("Transition Table: {}".format(TransitionTable))
    pl("Stranded Transition Table: {}".format(StrandedTTable))
    TransHandle.close()
    PileupHandle.close()
    return PileupTsv


def CustomPileupToTsv(inputBAM,
                      PileupTsv="default",
                      TransitionTable="default",
                      StrandedTTable="default",
                      bedfile="default",
                      progRepInterval=1000):
    '''
    A pileup tool for creating a tsv for each position in the bed file.
    Used for calling SNPs with high confidence.
    Also creates several tables:
    1. Counts for all consensus-->alt transitions (By Total and Merged reads)
    2. Counts for the above, specifying strandedness
    3. Number of Merged Reads supporting each allele
    '''
    TransTotalDict = {}
    TransMergedDict = {}
    StrandedTransTotalDict = {}
    StrandedTransMergedDict = {}
    NumTransitionsTotal = 0
    NumTransitionsMerged = 0
    MergedReadsProcessed = 0
    TotalReadsProcessed = 0
    import os.path
    if(PileupTsv == "default"):
        PileupTsv = inputBAM[0:-4] + '.Pileup.tsv'
    if(TransitionTable == "default"):
        TransitionTable = inputBAM[0:-4] + '.Trans.tsv'
    if(StrandedTTable == "default"):
        StrandedTTable = inputBAM[0:-4] + '.StrandedTrans.tsv'
    if(bedfile == "default"):
        return CustomPileupFullGenome(inputBAM, PileupTsv=PileupTsv,
                                      TransitionTable=TransitionTable,
                                      StrandedTTable="default",
                                      progRepInterval=progRepInterval)
    bedlines = HTSUtils.ParseBed(bedfile)
    if(os.path.isfile(inputBAM + ".bai") is False):
        pl("No bam index found for input bam - creating!")
        try:
            subprocess.check_call(['samtools', 'index', inputBAM])
        except subprocess.CalledProcessError:
            pl("Couldn't index BAM - coor sorting, then indexing!")
            inputBAM = HTSUtils.CoorSortAndIndexBam(inputBAM, uuid=True)
    NumPos = sum([line[2]-line[1] for line in bedlines])
    NumProcessed = 0  # Number of processed positions in pileup
    pl("Number of positions in bed file: {}".format(NumPos))
    bamHandle = pysam.AlignmentFile(inputBAM, "rb")
    PileupHandle = open(PileupTsv, "w")
    TransHandle = open(TransitionTable, "w")
    StrandedTransHandle = open(StrandedTTable, "w")
    FirstLine = True
    for line in bedlines:
        for pileupColumn in bamHandle.pileup(line[0], line[1], line[2]):
            NumProcessed += 1
            if((NumProcessed) % progRepInterval == 0):
                pl("Number of positions processed: {}".format(
                    NumProcessed + 1))
                pl("{:.1%} complete".format(NumProcessed / float(NumPos)))
                pl("Total reads processed: {}".format(TotalReadsProcessed))
                pl("Merged reads processed: {}".format(MergedReadsProcessed))
            PColSum = PCInfo(pileupColumn)
            MergedReadsProcessed += PColSum.MergedReads
            TotalReadsProcessed += PColSum.TotalReads
            if(FirstLine is True):
                PileupHandle.write(PColSum.toString(header=True))
                FirstLine = False
            else:
                PileupHandle.write(PColSum.toString())
            for key in PColSum.TransMergedCounts.keys():
                try:
                    TransMergedDict[
                        key] += PColSum.TransMergedCounts[key]
                except KeyError:
                    TransMergedDict[
                        key] = PColSum.TransMergedCounts[key]
                NumTransitionsMerged += PColSum.TransMergedCounts[
                    key]
            for key in PColSum.TransTotalCounts.keys():
                try:
                    TransTotalDict[
                        key] += PColSum.TransTotalCounts[key]
                except KeyError:
                    TransTotalDict[
                        key] = PColSum.TransTotalCounts[key]
                NumTransitionsTotal += PColSum.TransTotalCounts[
                    key]
            for key in PColSum.StrandedTransMergedCounts.keys():
                try:
                    StrandedTransMergedDict[
                        key] += PColSum.StrandedTransMergedCounts[
                            key]
                except KeyError:
                    StrandedTransMergedDict[
                        key] = PColSum.StrandedTransMergedCounts[
                            key]
            for key in PColSum.StrandedTransTotalCounts.keys():
                try:
                    StrandedTransTotalDict[
                        key] += PColSum.StrandedTransTotalCounts[
                            key]
                except KeyError:
                    StrandedTransTotalDict[
                        key] = PColSum.StrandedTransTotalCounts[
                            key]
    TransHandle.write(("Transition\tTotal Reads With Transition (Unflattened)"
                       "\tMerged Reads With Transition\tFraction Of Total "
                       "Transitions\tFraction Of Merged Transitions\n"))
    for key in TransTotalDict.keys():
        if(key[0] != key[3]):
            TransHandle.write("{}\t{}\t{}\n".format(key,
                              TransTotalDict[key],
                              TransMergedDict[key],
                              TransTotalDict[key]/float(
                                  NumTransitionsTotal),
                              TransMergedDict[key]/float(
                                  NumTransitionsMerged)))
    StrandedTransHandle.write(("Transition+Strandedness\tTotal Reads "
                               "(Unflattened)\tMergedReads With Transition\t"
                               "Fraction Of Total (Unflattened) Transitions"
                               "\tFraction of Merged Transitions\n"))
    for key in StrandedTransTotalDict.keys():
        if(key.split("&&")[0].split("->")[0] != key.split(
                "&&")[0].split("->")[1]):
            StrandedTransHandle.write("{}\t{}\t{}\n".format(key,
                                      StrandedTransTotalDict[key],
                                      StrandedTransMergedDict[key],
                                      StrandedTransTotalDict[key] / float(
                                          NumTransitionsTotal),
                                      StrandedTransMergedDict[key] / float(
                                          NumTransitionsMerged),
                                      ))
    pl("Transition Table: {}".format(TransitionTable))
    pl("Stranded Transition Table: {}".format(StrandedTTable))
    TransHandle.close()
    PileupHandle.close()
    return PileupTsv


def AlleleFrequenciesByBase(inputBAM, outputTsv="default",
                            progRepInterval=1000):
    if(outputTsv == "default"):
        outputTsv = inputBAM[0:-4] + '.allele.freq.tsv'
    if(os.path.isfile(inputBAM + ".bai") is False):
        pl("No bam index found for input bam - creating!")
        try:
            subprocess.check_call(['samtools', 'index', inputBAM])
        except subprocess.CalledProcessError:
            pl("Couldn't index BAM - coor sorting, then indexing!")
            inputBAM = HTSUtils.CoorSortAndIndexBam(inputBAM, uuid=True)
    NumProcessed = 0  # Number of processed positions in pileup
    inHandle = pysam.AlignmentFile(inputBAM, "rb")
    outHandle = open(outputTsv, "w")
    for pileup in inHandle.pileup():
        NumProcessed += 1
        if(NumProcessed % progRepInterval == 0):
            pl("Number of base positions processed: {}".format(NumProcessed))
        PColInfo = PCInfo(pileup)
    #TODO: Have it write contig, position, ref, and then counts for A, C, G, T
    # and then frequencies for A C G and T at each position, both by total
    # and merged counts.
    raise HTSUtils.ThisIsMadness("This script has not yet been completed.")
    return outputTsv

def MPileup(inputBAM, ref,
            bed="default",
            outputBCF="default",
            minbqual="20",
            minmqual="10"):
    if(outputBCF == "default"):
        if(len(inputBAM.split('.')) >= 6):
            outputBCF = inputBAM.split('.')[0] + ".fullMP.vcf"
        else:
            outputBCF = '.'.join(inputBAM.split('.')[0:-1]) + ".fullMP.vcf"
    if(bed != "default"):
        cmd = ("samtools mpileup -f {} -F 0.0000001 ".format(ref) +
               "-g -R -q " + minmqual + " -Q " + minbqual +
               " -l {} {}".format(bed, inputBAM) +
               " | bcftools view - > {}".format(outputBCF))
    else:
        cmd = ("samtools mpileup -f {} -F 0.0000001 ".format(ref) +
               "-I -S -g -D -R -q " + minmqual + " -Q " + minbqual +
               " " + inputBAM + " | bcftools view - > {}".format(outputBCF))
    pl("{} is command string".format(cmd))
    subprocess.check_call(cmd, shell=True)
    return outputBCF


def ParseVCF(inputVCFName):
    infile = open(inputVCFName, "r")
    VCFLines = [entry.strip().split('\t') for entry in infile.readlines(
    ) if entry[0] != "#" and entry.strip().split('\t')[4] != "<X>"]
    infile.seek(0)
    VCFHeader = [entry.strip(
    ) for entry in infile.readlines() if entry[0] == "#"]
    VCFEntries = [VCFRecord(
        line, inputVCFName) for line in VCFLines]
    ParsedVCF = VCFFile(VCFEntries, VCFHeader, inputVCFName)
    return ParsedVCF


def VCFRecordTest(inputVCFRec, filterOpt="default", param="default"):
    lst = [i.lower() for i in "bed,I16".split(',')]
    # print("lst = {}".format(lst))
    if(filterOpt.lower() not in lst):
        raise ValueError(("Filter option not supported. Available options: " +
                          ', '.join(lst)))
    passRecord = True
    if(filterOpt == "default"):
        raise ValueError("Filter option required.")
    if(filterOpt == "bed"):
        if(param == "default"):
            raise ValueError("Bedfile req. for bed filter.")
        bedReader = open(param, 'r')
        bedEntries = [l.strip().split('\t') for l in bedReader.readlines()]
        chr, pos = inputVCFRec.CHROM, inputVCFRec.POS
        chrMatches = [ent for ent in bedEntries if ent[0] == chr]
        try:
            posMatches = [match for match in chrMatches if match[
                          2] + 1 >= pos and match[1] + 1 <= pos]
            if len(posMatches) >= 1 and passRecord is True:
                passRecord = True
            else:
                passRecord = False
        except ValueError:
            raise ValueError("Malformed bedfile.")
            # return False
    # Set param to int, where it is the minimum dissent reads
    if(filterOpt == "I16"):
        if(param == "default"):
            raise ValueError("Minimum # dissenting reads must be set.")
        param = int(param)
        ConsensusIsRef = True
        I16Array = np.array(inputVCFRec.InfoArrayDict['I16']).astype("int")
        if(np.sum(I16Array[0:2]) < np.sum(I16Array[2:4])):
            ConsensusIsRef = False
        if(ConsensusIsRef is True):
            if(np.sum(I16Array[2:4]) > param):
                return True
            else:
                return False
        else:
            if(np.sum(I16Array[0:2]) > param):
                return True
            else:
                return False
    return passRecord
