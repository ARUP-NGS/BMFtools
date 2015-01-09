import argparse

from HTSUtils import printlog as pl
"""
bmftools contains various utilities for barcoded reads and for
somatic variant calling. Written to be in similar form to bcftools
and samtools.
"""

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="bmfsuites")
    SNVParser = subparsers.add_parser("snv")
    SNVParser.add_argument("inBAM",
                           help="Input BAM, Coordinate-sorted and indexed, "
                           "with BMF Tags included. bmftools runs on unflat"
                           "tened BAMs, but the results are not reliable be"
                           "cause this program assumes as much.")
    SNVParser.add_argument("--ref",
                           "-r",
                           help="Path to bwa index/ "
                           "faidx-indexed reference fasta.")

    args = parser.parse_args()
    return


if(__name__ == "__main__"):
    main()
