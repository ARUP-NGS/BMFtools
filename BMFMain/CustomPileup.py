#!/usr/bin/env python

import argparse

from BMFUtils.HTSUtils import printlog as pl


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--bam',
        '-i',
        help=("Coordinate-Sorted BAM File "),
        required=True
        )
    parser.add_argument(
        '--bed',
        '-b',
        help="Bedfile to use.",
        required=True)
    args = parser.parse_args()
    return


if(__name__ == "__main__"):
    main()
