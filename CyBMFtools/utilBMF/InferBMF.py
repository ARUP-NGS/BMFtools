#!/usr/bin/env python

import argparse

from utilBMF.HTSUtils import printlog as pl
from MawCluster.FauxBMF import *
from MawCluster.BCBam import ConsolidateInferred


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'inBAM',
        help="Input BAM file, name-sorted, with no supplementary "
        "or secondary alignments.")
    args = parser.parse_args()
    MarkedBAMDict = MarkReadPairPositions(args.inBAM)
    MarkedBAM = MarkedBAMDict["outBAM"]
    contigSets = MarkedBAMDict["contigSets"]
    RPSortBAM = SortBamByRPTag(MarkedBAM, contigSets=contigSets)
    ConsolidatedInferred = ConsolidateInferred(RPSortBAM)
    return


if(__name__ == "__main__"):
    main()
