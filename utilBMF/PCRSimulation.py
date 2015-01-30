#!/usr/bin/env
import argparse
"""
This program is motivated by the problem of deducing the effects of
PCR duplication and its error contributions in
Massively Parallel Sequencing datasets.
"""


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-n",
                        "--size-duplicate-family",
                        help="Number of duplications in the duplicate "
                        "family, one less than the number in the family",
                        default=10)
    parser.add_argument(
        "-i",
        "--iterations",
        help="Number of Monte Carlo iterations.",
        default=100000)
    parser.add_argument(
        "-r",
        "--random",
        help="Randomly select the number of reads for a duplicate family",
        default=False)
    parser.add_argument(
        "-o",
        "--outfile",
        help="File to write output.",
        default="")
    args = parser.parse_args()

    size = int(args.size_duplicate_family)
    outfile = args.outfile
    if(outfile == ""):
        outfile = "PCR_Duplication_Simulation_{}_InFamily_{}_Iterations.txt".format(
            size,
            args.iterations)
    countMC = 0
    sumHist = [0] * size
    print("Beginning loop")
    while(countMC < int(args.iterations)):
        iterationArray = dupDistribution(size)
        for count, entry in enumerate(iterationArray):
            sumHist[count] += entry
        countMC += 1
    outFile = open(outfile, "w", 0)
    outFile.write("@Description: PCR Duplication Simulation Results.\n")
    outFile.write(
        "@Header-Copied_#_times\tNumberOfReadsMatchingSaidDescription\n")
    for count, dup in enumerate(sumHist):
        outFile.write("{}\t{}\n".format(count, dup))
    outFile.write(repr(sumHist))
    return


def dupDistribution(size):
    import random
    random.seed()
    dupHist = [0.] * size
    # dupHist[0],dupHist[1]=1,1
    dupHist[0] = 1
    while(sum(dupHist) < size + 1):
        cdf = []
        for count, i in enumerate(dupHist):
            cdf.append(sum(dupHist[0:count]))
        sum_cdf = cdf[-1]
        bins = [(dupHist[x]) / sum_cdf for x in range(size)]
        dropTheNeedle = random.random()
        for count, bin in enumerate(bins):
            if(dropTheNeedle < bin):
                dupHist[count + 1] += 1
                break
    return dupHist

    print("Completing distribution of duplication.")
    print(repr(dupHist))
    return dupHist

if(__name__ == "__main__"):
    main()
