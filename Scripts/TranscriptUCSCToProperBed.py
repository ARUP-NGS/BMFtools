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
    parser.add_argument(
        '--bed',
        '-o',
        help="Output bed location. Default based on ucsc input.")
    args = parser.parse_args()
    header, entries = loadFile(args.ucsc)
    
    return


if(__name__ == "__main__"):
    main()


def loadFile(filename):
    infile = open(filename, "r")
    header = infile.readline().strip()
    firstPass = [line.strip().split('\t') for line in infile.readlines()]
    entries = [line for line in firstPass if 'NR' not in line[4]]
    return header, entries
