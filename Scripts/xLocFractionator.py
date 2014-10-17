import argparse
from collections import Counter

from HTSUtils import printlog as pl
from __builtin__ import enumerate


def main():
    pl("Now starting xLocFractionator")
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--genelist',
        '-g',
        help="Provide your file with a gene name on each line.")
    parser.add_argument(
        '--bed',
        '-b',
        help="Bed file with exon locations.",
        required=True)
    parser.add_argument('--outbed',
                        help="Output. Default based on genelist filename.",
                        default="default")
    parser.add_argument(
        "--ref",
        '-r',
        help="reference fasta",
        default="/mounts/genome/human_g1k_v37.fasta")
    parser.add_argument(
        '--fasta',
        '-o',
        help="output fasta",
        default="default")
    args = parser.parse_args()
    outbed = args.outbed
    outfasta = args.fasta
    if(outbed == "default"):
        outbed = '.'.join(args.genelist.split('.')[0:-1]) + '.cust.bed'
    if(outfasta == "default"):
        outfasta = '.'.join(args.bed.split('.')[0:-1]) + '.cust.fasta'
    gl = open(args.genelist, "r")
    genelist = [line.strip().lower() for line in gl]
    gl.close()
    pl(args.bed + " is args.bed")
    sourceBed = open(args.bed, "r")
    print("About to write to args.outbed: " + outbed)
    out = open(outbed, "w")
    bedData = [line.strip() for line in sourceBed.readlines()]
    print("Length of bedData is {}".format(len(bedData)))
    bedLines = [
        line +
        "\n" for line in bedData if line.split('\t')[4].lower() in genelist]
    startLines = ['\t'.join(line.split('\t')[0:3]) for line in bedLines]
    nameLines = ['.'.join(line.split('\t')[3:5]) for line in bedLines]
    endLines = ['\t'.join(line.split('\t')[6:]) for line in bedLines]
    unitedLines = []
    for count in range(len(startLines)):
        unitedLines.append(startLines[count] + "\t" + nameLines[
                           count] + "\t" + endLines[count])
    '''unitedBedLines = [line.split('\t')[0:2] +
                      [(line.split('\t')[3] +
                        '.' +
                        line.split('\t')[4])] +
                      line.split('\t')[5:] for line in bedLines]
    '''
    uniqTransDict = Counter([line.split('\t')[2] for line in unitedLines])
    renamedBedEntries = []
    print(uniqTransDict.keys())
    for key in uniqTransDict.keys():
        bedEntries = [l for l in unitedLines if l.split('\t')[2] == key]
        start = ['\t'.join(l.split('\t')[0:4]) for l in bedEntries]
        end = ['\t'.join(l.split('\t')[4:]) for l in bedEntries]
        for num, stLine in enumerate(start):
            stLine += '.Exon' + str(num+1)
            renamedBedEntries.append(stLine + '\t' + end[num])
    print(renamedBedEntries)
    out.writelines(renamedBedEntries)
    out.close()
    from subprocess import call
    call("wc -l {}".format(outbed), shell=True)
    pl("Now I am trying.")
    makeFasta(outbed, ref=args.ref, output=outfasta)
    return


def makeFasta(bedFile, ref="default", output="default"):
    print("Now grabbing from this bed: " + bedFile)
    from subprocess import call
    if(ref == "default"):
        raise ValueError("Required: fasta reference.")
    if(output == "default"):
        output = '.'.join(bedFile.split('.')[0:-1]) + '.custom.fasta'
    comStr = "fastaFromBed -name -fi {} -bed {} -fo".format(ref, bedFile)
    comStr += " " + output
    call(comStr, shell=True)
    return output


if(__name__ == "__main__"):
    main()
