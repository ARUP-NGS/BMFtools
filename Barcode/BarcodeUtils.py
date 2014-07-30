from Bio import SeqIO
import argparse
from pysam import Samfile

def generateHammingMatrix(input, outputTSV="default"):
    if(outputTSV=="default"):
        outputTSV = input.split('.')[0]+'.HammingMatrix.tsv'
    print("Output file for tab-separated values is {}.".format(outputTSV))
    import csv
    tsvWriter = open(outputTSV,'w',0)
    MatrixWriter = csv.writer(tsvWriter,delimiter="\t")
    inputBAM = Samfile(input,"rb")
    Headers = [entry.qname for entry in inputBAM]
    Headers.insert(0, "Names of the reads") 
    MatrixWriter.writerow(Headers)
    del Headers
    inputBAM = Samfile(input,"rb")#For iterating through again
    for entry in inputBAM:
        iteratorBAM = Samfile(input,"rb") #Opening multiple times to enable the looping
        BSentry = entry.tags[[i for i,j in enumerate(entry.tags) if j[0]=="BS"][0]][1]
        HammingDistanceRow = [hamming(BSentry, entry2.tags[[i for i,j in enumerate(entry2.tags) if j[0]=="BS"][0]][1]) for entry2 in iteratorBAM] 
        HammingDistanceRow.insert(0,entry.qname)
        MatrixWriter.writerow(HammingDistanceRow)    #^^A messy list comprehension within a list comprehension. It just returns the values of the hamming distances between the read "entry 1" and all of the reads in the BAM file
    inputBAM.close()
    tsvWriter.close()
    return outputTSV

def generateLevenshteinMatrix(input, outputTSV="default"):
    import Levenshtein
    if(outputTSV=="default"):
        outputTSV = input.split('.')[0]+'.HammingMatrix.tsv'
    print("Output file for tab-separated values is {}.".format(outputTSV))
    import csv
    tsvWriter = open(outputTSV,'w',0)
    MatrixWriter = csv.writer(tsvWriter,delimiter="\t")
    inputBAM = Samfile(input,"rb")
    Headers = [entry.qname for entry in inputBAM]
    Headers.insert(0, "Names of the reads") 
    MatrixWriter.writerow(Headers)
    del Headers
    inputBAM = Samfile(input,"rb")#For iterating through again
    for entry in inputBAM:
        iteratorBAM = Samfile(input,"rb") #Opening multiple times to enable the looping
        BSentry = entry.tags[[i for i,j in enumerate(entry.tags) if j[0]=="BS"][0]][1]
        HammingDistanceRow = [Levenshtein.distance(BSentry, entry2.tags[[i for i,j in enumerate(entry2.tags) if j[0]=="BS"][0]][1]) for entry2 in iteratorBAM] 
        HammingDistanceRow.insert(0,entry.qname)
        MatrixWriter.writerow(HammingDistanceRow)    #^^A messy list comprehension within a list comprehension. It just returns the values of the hamming distances between the read "entry 1" and all of the reads in the BAM file
    inputBAM.close()
    tsvWriter.close()
    return outputTSV

def hamming(str1, str2):
    import operator
    import Levenshtein
    from itertools import imap
    '''
    If you want to be flexible and allow different lengths, uncomment this block
    try:
        assert len(str1) == len(str2)
    except AssertionError:
        print("String 1 is {} with length {}, while String 2 is {} with length {}".format(str1,len(str1),str2,len(str2)))
        print("Calculating Levenshtein distance instead.")
        return Levenshtein.distance(str1,str2)
     '''   
    #ne = str.__ne__  ## this is surprisingly slow
    #ne = operator.ne
    return sum(imap(operator.ne, str1, str2))

#This function is to handle StopIterations with a little elegance
def has_elements(iterable):
    from itertools import tee
    iterable, any_check = tee(iterable)
    try:
        any_check.next()
        return True, iterable
    except StopIteration:
        return False, iterable