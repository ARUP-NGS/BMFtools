import argparse
import logging

from Bio import SeqIO
import pysam
from HTSUtils import printlog as pl


def generateHammingMatrix(input, outTSV="default"):
    if(outTSV == "default"):
        outTSV = input.split('.')[0] + '.HammingMatrix.tsv'
    print('Output file for tab-separated values is {}.'.format(outTSV))
    import csv
    tsvWriter = open(outTSV, 'w', 0)
    MatrixWriter = csv.writer(tsvWriter, delimiter="\t")
    inputBAM = pysam.Samfile(input, "rb")
    Headers = [entry.qname for entry in inputBAM]
    Headers.insert(0, "Names of the reads")
    MatrixWriter.writerow(Headers)
    del Headers
    inputBAM = pysam.Samfile(input, "rb")  # For iterating through again
    for entry in inputBAM:
        iteratorBAM = pysam.Samfile(input, "rb")  # Opening multiple times
        BSentry = entry.tags[[i for i,
                              j in enumerate(
                                  entry.tags) if j[0] == "BS"][0]][1]
        HDisRow = [hamming(BSentry, entry2.tags[[i for i,
                   j in enumerate(entry2.tags) if j[
                       0] == "BS"][0]][
                   1]) for entry2 in iteratorBAM]
        HDisRow.insert(0, entry.qname)
        MatrixWriter.writerow(HDisRow)
    inputBAM.close()
    tsvWriter.close()
    return outTSV

'''
def generateLevenshteinMatrix(input, outTSV="default"):
    import Levenshtein
    if(outTSV == "default"):
        outTSV = input.split('.')[0] + '.HammingMatrix.tsv'
    pl("Output file for tab-separated values is {}.".format(outTSV))
    import csv
    tsvWriter = open(outTSV, 'w', 0)
    MatrixWriter = csv.writer(tsvWriter, delimiter="\t")
    inputBAM = pysam.Samfile(input, "rb")
    Headers = [entry.qname for entry in inputBAM]
    Headers.insert(0, "Names of the reads")
    MatrixWriter.writerow(Headers)
    del Headers
    inputBAM = pysam.Samfile(input, "rb")  # For iterating through again
    for entry in inputBAM:
        iteratorBAM = pysam.Samfile(input, "rb")  # Opening each time
        BSentry = entry.tags[[i for i,
                              j in enumerate(entry.tags) if j[
                                  0] == "BS"][0]][1]
        HDisRow = [Levenshtein.distance(BSentry, entry2.tags[
            [i for i, j in enumerate(
                entry2.tags) if j[0] == "BS"][
                0]][1]) for entry2 in iteratorBAM]
        HDisRow.insert(0, entry.qname)
        MatrixWriter.writerow(HDisRow)
    inputBAM.close()
    tsvWriter.close()
    return outTSV
'''


def hamming(str1, str2):
    import operator
    from itertools import imap
    '''
    If you want to allow different lengths, uncomment this block
    import Levenshtein
    try:
        assert len(str1) == len(str2)
    except AssertionError:
        str="String 1 is {} with length {}".format(str1,len(str1))
        str+=" while String 2 iswith length {}".format(str2, len(str2))
        pl(str)
        pl("Calculating Levenshtein distance instead.")
        return Levenshtein.distance(str1,str2)
     '''
    # ne = str.__ne__  ## this is surprisingly slow
    # ne = operator.ne
    return sum(imap(operator.ne, str1, str2))


# This function is to handle StopIterations with a little elegance
def has_elements(iterable):
    from itertools import tee
    iterable, any_check = tee(iterable)
    try:
        any_check.next()
        return True, iterable
    except StopIteration:
        return False, iterable
