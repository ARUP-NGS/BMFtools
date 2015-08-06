"""
Contains functions for analyzing insert size in the context of variant allele
support.
"""

import cython
import numpy as np
import pysam
from utilBMF.HTSUtils import pPileupRead
from utilBMF.ErrorHandling import ThisIsMadness as Tim, ImproperArgumentError
from utilBMF.Inliners import Num2NucN
from array import array as py_array


cdef class SNVAlleleWrangler:


    cdef inline bint pass_record(self, pPileupRead_t PR, int nuc_index):
        if Num2NucN(nuc_index) != <char>PR.BaseCall:
            return False
        if (PR.FM < self.filter_cnf[0] or PR.MQ < self.filter_cnf[1]
            or PR.BQ < self.filter_cnf[2] or PR.FA < self.filter_cnf[3]
            or <double_t> PR.FA / PR.FM < self.filter_cnf[4]):
            return False
        return True


    cdef inline void populate(self):
        cdef pPileupRead_t PR
        cdef int nuc_index
        self.alleles = {nuc_index: [PR for PR in self.column if
                                    self.pass_record(PR, nuc_index)]
                        for nuc_index in xrange(5)}


    cdef inline void fast_forward(self):
        """Move forward until we arrive at the correct position"""
        cdef IteratorColumnRegion_t iterator
        cdef PileupColumn_t column
        iterator = self.handle.pileup(self.handle.getrname(self.reference_id),
                                      self.pos)
        column = iterator.next()
        try:
            assert column.pos <= self.pos
            while column.pos < self.pos:
                column = iterator.next()
        except StopIteration:
            raise Tim("Position not in BAM - are we sure that this "
                      "position is even on this contig?")
        self.column = column

    def __init__(self, cystr bampath, cystr contig, int32_t pos,
                 int32_t minFM=0, int32_t minMQ=0, int32_t minPV=0,
                 int32_t minFA=0, double_t minFAFrac=0.0):
        # C definitions
        cdef dict tmpdict

        self.handle = pysam.AlignmentFile(bampath, "rb")
        self.pos = pos
        self.reference_id = self.handle.gettid(contig)

        self.filter_cnf = py_array('f', [<double_t>minFM, <double_t>minMQ,
                                         <double_t> minPV, <double_t> minFA,
                                         minFAFrac])

        # Test to see if this position is even available on the contig.
        '''
        try:
            tmpdict = self.handle.header
            assert tmpdict['SQ'][self.reference_id]['LN'] > pos
        except AssertionError:
            raise ImproperArgumentError("Position selected is greater than"
                                        " the length of contig. Abort!")
        '''
        self.fast_forward()
        self.populate()


