import argparse

from HTSUtils import printlog as pl


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--ucsc',
        '-i',
        help=(
            "Provide your file 2nd and 3rd columns containing "
            "comma-delimited exon starts and stops "
            "respectively.")
        )
    args = parser.parse_args()
    return


if(__name__ == "__main__"):
    main()
