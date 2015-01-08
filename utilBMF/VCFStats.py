#!/usr/bin/env python
import argparse
import logging

from MawCluster.BCVCF import VCFStats
from HTSUtils import printlog as pl

# Global Variables
Logger = logging.getLogger("Primarylogger")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'inVCF',
        help=("Coordinate-sorted, indexed BAM"),
        )
    args = parser.parse_args()
    print("Beginning VCFStats")
    VCFStatsTable = VCFStats(args.inVCF)
    print("Stats table path: {}".format(VCFStatsTable))
    return VCFStatsTable

if(__name__ == "__main__"):
    main()
