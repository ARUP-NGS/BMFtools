#!/mounts/anaconda/bin/python

import argparse

from MawCluster import BCVCF

#  Cleans up FreeBayes/GATK VCFEntries which should be split into
#  multiple lines.


def fixVCF(inputVCF, outputVCF="default", write="False"):
    if(outputVCF == "default"):
        outputVCF = '.'.join(inputVCF.split('.')[0:-1]) + '.fixed.vcf'
    startVCF = BCVCF.ParseVCF(inputVCF)
    outVCFHeader = startVCF.header
    stdEntries = [i for i in startVCF.Records if(
                  len(i.REF) != len(i.ALT) or len(i.REF) == 1)]
    elseEntries = [i for i in startVCF.Records if(
                   len(i.REF) == len(i.ALT) and len(i.REF) > 1)]
    NewEntries = []
    for e in elseEntries:
        for num, pair in enumerate(zip(e.REF, e.ALT)):
            e.REF, e.ALT = pair[0], pair[1]
            if(num != 0):
                e.POS = str(int(e.POS) + 1)  # Incrementing location
            newRecord = BCVCF.VCFRecord(e.toString().split('\t'), outputVCF)
            if(newRecord.ALT != newRecord.REF):
                stdEntries.append(newRecord)
    for s in stdEntries:
        NewEntries.append(s)
    outVCF = BCVCF.VCFFile(NewEntries, outVCFHeader, outputVCF)
    if(write.lower() == "true"):
        outVCF.write(outputVCF)
    return outputVCF, outVCF


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i',
        '--vcf',
        help="Provide your vcf file.",
        required=True)
    parser.add_argument(
        '-f',
        '--force',
        help="Set to True to replace original VCF file.",
        default="False")
    parser.add_argument(
        '-o',
        '--output',
        help="Provide filename/path for output file. Optional")
    args = parser.parse_args()
    print("Now beginning cleaning of file {}".format(args.vcf))
    outVCFFile, outVCFObj = fixVCF(args.vcf)
    if(args.force.lower() == "true"):
        print("Overwriting original file.")
        outVCFObj.write(args.vcf)
    else:
        outVCFObj.write(outVCFFile)
    return

if(__name__ == "__main__"):
    main()
