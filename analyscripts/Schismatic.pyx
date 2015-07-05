# cython: cdivision=True
# Standard library imports
from itertools import groupby
from array import array

# Third party imports
import numconv
import pysam

# BMFTools imports
from utilBMF.HTSUtils import pFastqProxy, pFastqFile, getBS, RevCmp


converter = numconv.NumConv(4, "ACGT")
int2nuc = [converter.int2str(x) for x in xrange(5000)]


def SplitFastqFamilies(cystr inFq, cystr outfqpath,
                       int FamilyCap):
    """
    Splits up fastq families into subfamilies of size
    FamilyCap by appending nucleotides to it.
    """
    cdef pFastqFile_t inHandle = pFastqFile(inFq)
    cdef pFastqProxy_t read
    cdef cystr barcode
    cdef object generator, ohw
    cdef size_t MemberCount
    cdef size_t SubFamilyCount
    outHandle = open(outfqpath, "w")
    ohw = outHandle.write
    for barcode, generator in groupby(inHandle, GetFqBS):
        SubFamilyCount = 0
        for MemberCount, read in enumerate(generator):
            read.comment = read.comment.replace(
                "BS=" + barcode, "BS=" + barcode + int2nuc[SubFamilyCount])
            ohw(str(read))
            if(MemberCount % FamilyCap == 0):
                SubFamilyCount += 1
    return outfqpath
