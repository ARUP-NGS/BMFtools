#/mounts/anaconda/bin/python

from BarcodeHTSTools import IllegalArgumentError

#TODO: Write VCFRecordTest
#TODO: 

class VCFFile:
    """A simple VCFFile object, consisting of a header, a name for the file from which they came, and a list of all VCFRecords."""
    def __init__(self,VCFEntries,VCFHeader,inputVCFName):
        self.sampleName = inputVCFName
        self.header = VCFHeader
        self.Records = VCFEntries
        self.sampleNamesArray = [inputVCFName]
        self.numSamples = len(self.sampleNamesArray)
        
    def cleanRecords(self):
        NewRecordList = [entry for entry in self.Records if entry.ALT != "X"]
        self.Records = NewRecordList
       
    def filter(self,filterOpt="default",bed="default"):
        if(filterOpt=="default"):
            try:
                raise IllegalArgumentError("Sorry, I can't filter your VCF if you don't provide a filter method")
            except IllegalArgumentError:
                print("Returning nothing.")
                return
        NewVCFEntries = []
        for entry in self.VCFEntries:
                if(VCFRecordTest(entry, filterOpt)==True):
                    NewVCFEntries.append(entry)
        if(filterOpt=="bed"):
            if(bed=="default"):
                try:
                    raise IllegalArgumentError("I can't filter by a bed file if you don't provide one! (Fine, I'll take a Large, but I don't GAF, alright??")
                except IllegalArgumentError:
                    print("Returning nothing.")
            filterOpt = filterOpt + "&" +  bed
        NewVCFFile = VCFFile(NewVCFEntries,self.header,self.sampleName + "FilteredBy{}".format(filterOpt))
        #TODO: make a new VCFFile object based on location.
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
            raise AttributeError("This VCFFile object has not been initialized.")
        
    def write(self,Filename):
        FileHandle = open(Filename,"w")
        for headerLine in self.header:
            FileHandle.write(headerLine)
        for VCFEntry in self.Records:
            recordStr = '\t'.join([VCFEntry.CHROM,VCFEntry.POS,VCFEntry.ID,VCFEntry.REF,VCFEntry.ALT,VCFEntry.QUAL, \
                                 VCFEntry.FILTER,VCFEntry.INFO,VCFEntry.FORMAT,VCFEntry.GENOTYPE,'\t'.join(VCFEntry.Samples)])
            FileHandle.write(recordStr)
        FileHandle.close()
        return

class VCFRecord:
    """A simple VCFRecord object, taken from an item from the list which ParseVCF returns as "VCFEntries" """
    def __init__(self,VCFEntry,VCFFilename):
        self.CHROM = VCFEntry[0]
        self.POS = VCFEntry[1] 
        self.ID = VCFEntry[2]
        self.REF = VCFEntry[3]
        self.ALT = VCFEntry[4] 
        self.QUAL = VCFEntry[5]
        self.FILTER = VCFEntry[6]
        self.INFO = VCFEntry[7]
        self.InfoKeys=[entry.split('=')[0] for entry in self.INFO.split(';')]
        self.InfoValues=[entry.split('=')[1] for entry in self.INFO.split(';')]
        self.InfoValArrays = [entry.split(',') for entry in self.InfoValues]
        self.InfoDict=dict(zip(self.InfoKeys,self.InfoValues)) 
        self.InfoArrayDict=dict(zip(self.InfoKeys,self.InfoValArrays)) 
        self.FORMAT = VCFEntry[8]
        self.GENOTYPE = VCFEntry[9]
        self.GenotypeDict = dict(zip(self.FORMAT.split(':'),self.GENOTYPE.split(':')))
        self.GenotypeKeys = self.FORMAT.split(':')
        self.GenotypeValues = self.GENOTYPE.split(':')
        self.Samples = [""]
        if(len(VCFEntry) > 10):
            for field in VCFEntry[10:]:
                self.Samples.append(field)
        self.VCFFilename = VCFFilename
        
    def update(self):
        self.InfoValues = [','.join(InfoValArray) for InfoValArray in self.InfoValArrays]
        infoEntryArray = [InfoKey+"="+InfoValue for InfoKey, InfoValue in zip(self.InfoKeys, self.InfoValues)]
        self.INFO = ';'.join(infoEntryArray)
        self.InfoKeys=[entry.split('=')[0] for entry in self.INFO.split(';')]
        self.InfoValues=[entry.split('=')[1] for entry in self.INFO.split(';')]
        self.InfoValArrays = [entry.split(',') for entry in self.InfoValues]
        self.InfoDict=dict(zip(self.InfoKeys,self.InfoValues)) 
        self.InfoArrayDict=dict(zip(self.InfoKeys,self.InfoValArrays)) 
        self.GenotypeKeys = self.FORMAT.split(':')
        self.GenotypeValues = self.GENOTYPE.split(':')
        self.FORMAT = ":".join(self.GenotypeKeys)
        self.GENOTYPE = ":".join(self.GenotypeValues)
        self.GenotypeDict = dict(zip(self.FORMAT.split(':'),self.GENOTYPE.split(':')))
            
    def str(self):
        self.update()
        if(len(self.Samples==0)):
            recordStr = '\t'.join([self.CHROM,self.POS,self.ID,self.REF,self.ALT,self.QUAL, \
                                 self.FILTER,self.INFO,self.FORMAT,self.GENOTYPE])
        return recordStr
        
#I also want to be able to grab all of the records for a given record, as well as grab the file from which the records came.

def CleanupPileup(inputPileup,outputPileup="default"):
    import subprocess
    if(outputPileup=="default"):
        outputPileup='.'.join(inputPileup.split('.')[0:-1]) + ".xrm.vcf"
    commandStr = "awk '$5!=\"X\"' {} > {}".format(inputPileup,outputPileup)
    subprocess.call(commandStr,shell=True)
    return outputPileup

def criteriaTest(entry,filterOpt="default",bed="default"):
    list="bed"
    if(filterOpt=="default"):
        print("List of valid filters: {}".format(', '.join(list)))
        raise ValueError("Filter must be set! Requires an exact match (case insensitive).")
    if(filterOpt=="bed"):
        if(bed=="default"):
            raise ValueError("Bed file must be provided in order to filter thereby!")
        bedHandle = open(bed,"r")
        bedEntries = [line.strip().split('\t') for line in bedHandle]
        for line in bedEntries:
            if(line[0]!=entry.CHROM):
                continue
            else:
                if(int(entry.POS)-1 >= int(line[1]) and int(entry.POS)-1 <= int(line[2])):
                    bedHandle.close()
                    return True
                else:
                    continue
        bedHandle.close()
        return False
    raise RuntimeWarning("This should never happen. It seems that no valid filters were set and then something went horribly wrong!")
    return False



def MPileup(inputBAM,bedfile,ref, outputBCF="default"):
    import subprocess
    if(outputBCF=="default"):
        if(len(inputBAM.split('.')) >= 6):
            outputBCF = inputBAM.split('.')[0] + ".fullMP.vcf";
        else:
            outputBCF = '.'.join(inputBAM.split('.')[0:-1]) + ".fullMP.vcf";
    commandStr = "samtools mpileup -f {} -F 0.0001 -I -S -g -D -R -q 10 -Q 30 -l {} {} | bcftools view - > {}".format(ref,bedfile,inputBAM,outputBCF)
    print("{} is command string".format(commandStr))
    subprocess.call(commandStr, shell=True)
    return outputBCF

def ParseVCF(inputVCFName):
    infile = open(inputVCFName,"r")
    VCFLines = [entry.strip().split('\t') for entry in infile.readlines() if entry[0]!="#"]
    VCFHeader = [entry.strip() for entry in infile.readlines() if entry[0]=="#"]
    VCFEntries = [VCFRecord(entry,inputVCFName) for entry in VCFLines]
    ParsedVCF = VCFFile(VCFEntries,VCFHeader,inputVCFName)
    return ParsedVCF
