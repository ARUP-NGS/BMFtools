#!/usr/bin/env python
import argparse

import cython
cimport cython
import numpy as np
cimport numpy as np
DTYPE = np.ndarray
ctypedef np.ndarray DTYPE_t

"""
This program is motivated by the problem of deducing the effects of
PCR duplication and its error contributions in
Massively Parallel Sequencing datasets.
"""


@cython.locals(countMC=cython.long, size=cython.long)
@cython.boundscheck(False)
@cython.wraparound(False)
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
        default=10)
    parser.add_argument(
        "-o",
        "--outfile",
        help="File to write output.",
        default="default")
    args = parser.parse_args()

    size = int(args.size_duplicate_family)
    outfile = args.outfile
    if(outfile == "default"):
        outfile = "PCR_Dup_Simulation_{}_InFamily_{}_Iterations.txt".format(
            size,
            args.iterations)
    countMC = 0
    cdef np.ndarray[np.int, ndim = 1] sumHist = np.zeros(size)
    print("Beginning analysis.")
    print("Number of copies: {}".format(size))
    print("Number of times to simulate: {}".format(args.iterations))
    while(countMC < int(args.iterations)):
        print("Beginning iteration {}".format(countMC + 1))
        sumHist = np.sum(sumHist, dupDistribution(size))
        countMC += 1
    outFile = open(outfile, "w", 0)
    outFile.write("#Description: PCR Duplication Simulation Results.\n")
    outFile.write(
        "#Times Copied\tNumberOfReadsMatchingSaidDescription\n")
    for count, dup in enumerate(sumHist):
        outFile.write("{}\t{}\n".format(count, dup))
    outFile.write(repr(sumHist))
    return


@cython.locals(size=cython.long, i=cython.long)
@cython.boundscheck(False)  # Turn off boundscheck for this function
@cython.wraparound(False)  # Turn off negative indexing
@cython.returns(DTYPE_t)
def dupDistribution(size):
    print("Beginning simulation for n={}".format(size))
    cdef np.ndarray[np.int, ndim = 1] dupHist = np.zeros(size)
    # dupHist[0],dupHist[1]=1,1
    dupHist[0] = 1
    while(sum(dupHist) < size + 1):
        '''
        for count, i in enumerate(dupHist):
            cdf.append(np.sum(dupHist[0:count]))
        '''
        bins = [np.sum(dupHist[0:i]) for i in range(size)] / np.max()
        try:
            dupHist[np.nonzero(bins > np.random.random())[0][0] + 1] += 1
            # Increases the number of molecules duplicated as many times as
            # that bin plus one.
        except IndexError:
            print(repr(dupHist))
            print(repr(bins))
            return 1
    print("Completing distribution of duplication.")
    print(repr(dupHist))
    return dupHist

if(__name__ == "__main__"):
    main()
