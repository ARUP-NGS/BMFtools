#!/mounts/anaconda/bin/python

class VCFFile:
    """A simple VCFFile object, consisting of a header, a name for the file from which they came, and a list of all VCFRecords."""
    def __init__(self,VCFEntries,VCFHeader,inputVCFName):
        self.sampleName = inputVCFName
        self.header = VCFHeader
        self.Records = VCFEntries
       
    def bedFilter(self,BedPath):
        #TODO: make a new VCFFile object based on location.
        return
    
    def len(self):
        try:
            return len(self.Records)
        except AttributeError:
            raise AttributeError("This VCFFile object has not been initialized.")
        
    def write(self,FileHandle):
        for headerLine in self.header:
            FileHandle.write(headerLine)
        #TODO: Finish this write function!

class VCFRecord:
    """A simple VCFRecord object, taken from an item from the list which ParseVCF returns as "VCFEntries" """
    def __init__(self,VCFEntry,VCFFilename):
        self.CHROM, = VCFEntry[0]
        self.POS = VCFEntry[1] 
        self.ID = VCFEntry[2]
        self.REF = VCFEntry[3]
        alt = VCFEntry[4].split(',')
        try:
            self.ALT = alt.remove("X")
        except ValueError:
            self.ALT = alt 
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
        self.Samples = []
        if(len(VCFEntry) > 10):
            for field in VCFEntry[10:]:
                self.Samples.append(field)
        self.VCFFilename = VCFFilename
                
    def view(self):
        return #TODO: finish this view/display function!
        
#I also want to be able to grab all of the records for a given record, as well as grab the file from which the records came.

def CleanupPileup(inputPileup,outputPileup="default"):
    import subprocess
    if(outputPileup=="default"):
        outputPileup='.'.join(inputPileup.split('.')[0:-1]) + ".xrm.vcf"
    commandStr = "awk '$5!=\"X\" {} > {}".format(inputPileup,outputPileup)
    subprocess.call(commandStr,shell=True)
    subprocess.call("bcftools index {}".format(outputPileup))
    return outputPileup

def MPileup(inputBAM,bedfile,ref, outputBCF="default"):
    import subprocess
    if(outputBCF=="default"):
        if(len(inputBAM.split('.')) >= 6):
            outputBCF = inputBAM.split('.')[0] + ".fullMP.vcf";
        else:
            outputBCF = '.'.join(inputBAM.split('.')[0:-1]) + ".fullMP.vcf";
    commandStr = "samtools mpileup -f {} -F 0.0001 -I -S -g -D -R -q 10 -Q 30 -l {} {} | bcftools view > {}".format(ref,bedfile,inputBAM,outputBCF)
    subprocess.call(commandStr, shell=True)
    return outputBCF

def ParseVCF(inputVCFName):
    infile = open(inputVCFName,"r")
    VCFEntries = [entry.strip().split('\t') for entry in infile.readlines() if entry[0]!="#"]
    VCFHeader = [entry.strip() for entry in infile.readlines() if entry[0]=="#"]
    ParsedVCF = VCFFile(VCFEntries,VCFHeader,inputVCFName)
    return ParsedVCF