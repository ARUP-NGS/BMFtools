import subprocess
import decimal
import numpy as np

from utilBMF.HTSUtils import ThisIsMadness, printlog as pl
from utilBMF import HTSUtils

"""
Contains tools for working with VCF Files - writing, reading, processing.
"""


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
        try:
            self.POS = VCFEntry[1]
        except IndexError:
            print("Line: {}".format(VCFEntry))
            raise ThisIsMadness("Something went wrong.")
        self.ID = VCFEntry[2]
        self.REF = VCFEntry[3]
        if("<X>" in VCFEntry[4]):
            self.ALT = VCFEntry[4].replace(",<X>", "")
        else:
            self.ALT = VCFEntry[4]
        """
        if("<X>" != VCFEntry[4]):
            self.ALT = ','.join(VCFEntry[4].split(',').remove("<X>"))
        else:
            self.ALT = "<X>"
        """
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
                self.InfoUnpaired.append(entry)
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
        except KeyError:
            pass
            # print("No I16 present; continuing.")
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


def is_reverse_to_str(boolean):
    if(boolean is True):
        return "reverse"
    elif(boolean is False):
        return "forward"
    else:
        return "unmapped"


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
    try:
        infile = open(inputVCFName, "r")
    except TypeError:
        raise TypeError("Argument provided: {}".format(inputVCFName))
    VCFLines = [entry.strip().split('\t') for entry in infile.readlines(
    ) if entry[0] != "#"]
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


def VCFStats(inVCF, TransCountsTable="default"):
    print("About to run VCFStats on {}".format(inVCF))
    if(TransCountsTable == "default"):
        TransCountsTable = inVCF[0:-3] + "trans.vcf.tsv"
    inVCF = ParseVCF(inVCF)
    TransCountsTableHandle = open(TransCountsTable, "w")
    TransitionDict = {}
    TransitionPASSDict = {}
    TransitionCountsDict = {}
    TransitionCountsPASSDict = {}
    RefConsCallsCountsDict = {}
    RefConsCallsCountsPASSDict = {}
    TransCountsTableHandle.write("Transition\tCount\tCount marked \"PASS"
                                 "\"\tFraction Of Calls With Given Ref/Cons\t"
                                 "MeanAlleleFraction\tFraction Of PASS Calls "
                                 "With Given Ref/Cons\tMean Allele Fraction "
                                 "Of PASS Calls\n")
    for RefCons in ["A", "C", "G", "T"]:
        for Var in ["A", "C", "G", "T"]:
            if(RefCons == Var):
                continue
            TransitionDict[RefCons + "-->" + Var] = [
                rec for rec
                in inVCF.Records
                if rec.InfoArrayDict['CONS'] == RefCons
                or rec.REF == RefCons and rec.ALT == Var]
            TransitionPASSDict[RefCons + "-->" + Var] = [
                rec for rec
                in inVCF.Records
                if (rec.InfoArrayDict['CONS'] == RefCons
                    or rec.REF == RefCons) is True and (rec.FILTER == "PASS"
                                                        and rec.ALT == Var)]
            TransitionCountsDict[RefCons + "-->" + Var] = len(
                TransitionDict[RefCons + "-->" + Var])
            TransitionCountsPASSDict[RefCons + "-->" + Var] = len(
                TransitionPASSDict[RefCons + "-->" + Var])
    for RefCons in ["A", "C", "G", "T"]:
        for key in TransitionCountsDict.keys():
            if key[0] == RefCons:
                try:
                    RefConsCallsCountsDict[
                        RefCons] += TransitionCountsDict[key]
                except KeyError:
                    RefConsCallsCountsDict[RefCons] = TransitionCountsDict[key]
                try:
                    RefConsCallsCountsPASSDict[
                        RefCons] += TransitionCountsPASSDict[key]
                except KeyError:
                    RefConsCallsCountsPASSDict[
                        RefCons] = TransitionCountsPASSDict[key]
    MeanAlleleFractionDict = {}
    MeanAlleleFractionPASSDict = {}
    for key in TransitionCountsDict.keys():
        MeanAlleleFractionDict[key] = np.mean([float(rec.InfoDict['AF']) for
                                               rec in TransitionDict[key] if
                                               float(
                                                   rec.InfoDict['AF']) < 0.1])
        MeanAlleleFractionPASSDict[key] = np.mean(
            [float(rec.InfoDict['AF'])for rec in
             TransitionPASSDict[key] if float(rec.InfoDict['AF']) < 0.1])
    TransitionFractionForRefConsDict = {}
    TransitionFractionForRefConsPASSDict = {}
    for key in TransitionCountsDict.keys():
        try:
            TransitionFractionForRefConsDict[
                key] = TransitionCountsDict[key] / float(
                    RefConsCallsCountsDict[key[0]])
        except ZeroDivisionError:
            TransitionFractionForRefConsDict[key] = 0
        try:
            TransitionFractionForRefConsPASSDict[
                key] = TransitionCountsPASSDict[key] / float(
                    RefConsCallsCountsPASSDict[key[0]])
        except ZeroDivisionError:
            TransitionFractionForRefConsPASSDict[key] = 0
    for key in TransitionCountsDict.keys():
        TransCountsTableHandle.write("\t".join(
            [str(i)[0:8] for i in [key,
                                   TransitionCountsDict[key],
                                   TransitionCountsPASSDict[key],
                                   TransitionFractionForRefConsDict[key],
                                   MeanAlleleFractionDict[key],
                                   TransitionFractionForRefConsPASSDict[key],
                                   MeanAlleleFractionPASSDict[key]
                                   ]]) + "\n")
    TransCountsTableHandle.close()
    return TransCountsTable
