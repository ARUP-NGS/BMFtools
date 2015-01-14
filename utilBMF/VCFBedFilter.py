#!/usr/bin/env python
import argparse
import logging

from MawCluster.BCVCF import FilterVCFFileByBed
from utilBMF.HTSUtils import printlog as pl

# Global Variables
Logger = logging.getLogger("Primarylogger")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'inVCF',
        help=("Input VCF"),
        )
    parser.add_argument(
        '--outVCF',
        "-o",
        help=("Output VCF"),
        default="default"
        )
    parser.add_argument(
        '--bed',
        "-b",
        help=("Bed file"),
        )
    args = parser.parse_args()
    inVCF = args.inVCF
    if(args.outVCF != "default"):
        outVCF = args.outVCF
    else:
        outVCF = inVCF[0:-4] + ".bedfilter.vcf"
    pl("Filtering by bed: {}".format(args.bed))
    pl("Output VCF: {}".format(outVCF))
    print("Beginning to filter VCF by bed file.")
    outVCF = FilterVCFFileByBed(inVCF, args.bed, outVCF=outVCF)
    print("Output VCF: {}".format(outVCF))
    return outVCF

if(__name__ == "__main__"):
    main()
