from BCFastq import halveFqRecords

import argparse

from HTSUtils import printlog as pl


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--fq',
        '-i',
        help=(
            "Provide your file 2nd and 3rd columns containing "
            "comma-delimited exon starts and stops "
            "respectively.")
    )
    parser.add_argument(
        '--outfq1',
        help="output fastq file",
        default="default")
    parser.add_argument(
        '--outfq2',
        help="output fastq file",
        default="default")
    parser.add_argument(
        '--minLength',
        '-m',
        help="Minimum length for a read to be halved and kept.")
    pl("Beginning halving of the fastq reads provided.")
    args = parser.parse_args()
    halveFqRecords(args.fq, outfq1=args.outfq1, outfq2=args.outfq2, minLength=40)
    return


if(__name__ == "__main__"):
    main()
