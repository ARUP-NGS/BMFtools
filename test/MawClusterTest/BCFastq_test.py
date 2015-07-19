#!/usr/bin/env python
import cython
import shlex
import subprocess
import sys
import unittest

from BMFMain.Workflow import pairedFastqShades
from utilBMF.HTSUtils import pFastqFile
from MawCluster.BCFastq import (pCompareFqRecsFast as cCompareFqRecsFast,
                                pairedFastqConsolidate, singleFastqConsolidate)

__author__ = 'dnephi and BrettKennedy'


class MyTestCase(unittest.TestCase):
    """
    TODO: Add a unit test to assert that the sum of the sizes of merged
    families equals the number of pre-merged reads.
    """

    def tearDown(self):
        self.handle1.close()
        self.handle2.close()
        for filename in self.filenames:
            subprocess.check_call(["rm", filename])

    @cython.locals(i=pFastqProxy_t)
    def setUp(self):
        self.handle1 = pFastqFile("TinyTest_R1.BS.fastq.gz")
        self.handle2 = pFastqFile("TinyTest_R2.BS.fastq.gz")
        self.prefastq1 = [i for i in self.handle1]
        self.prefastq2 = [i for i in self.handle2]
        self.filenames = []

    @cython.locals(PASS=bint)
    def test_FqProc(self):
        PASS = True

        if not PASS:
            raise AssertionError(
                "Str: '%s' does not match unit test. Abort!" % tmpStr)

    def test_pfc(self):
        pairedFastqShades("../data/TinyDemo_R1.fastq.gz",
                          "../data/TinyDemo_R3.fastq.gz",
                          indexFq="../data/TinyDemo_R2.fastq.gz")
        cStrToCheck = ("cat TinyDemo_R1.shaded.BS.cons.fastq | paste "
                       "- - - - | cut -d'|' -f4 | cut -d'=' -f2 | paste -sd+ "
                       "| bc")
        assert int(subprocess.check_output(cStrToCheck,
                                           shell=True).strip()) == 25000
        self.filenames += ["TinyDemo_R1.fastq.famstats.txt",
                           "TinyDemo_R1.shaded.BS.cons.fastq",
                           "TinyDemo_R1.shaded.BS.fastq",
                           "TinyDemo_R1.shaded.fastq",
                           "TinyDemo_R3.shaded.BS.cons.fastq",
                           "TinyDemo_R3.shaded.BS.fastq",
                           "TinyDemo_R3.shaded.fastq"]

    def test_sfc(self):
        singleFastqConsolidate("../data/TestR1.fastq")
        conFq = open("TestR1.cons.fastq")
        assert conFq.readlines() == ['@AAAAAATGGACCCATTAACC 1:N:0:ACAGTG|FP=I'
                                     'ndexPass|BS=AAAAAATGGACCCATTAACC|FM=3|N'
                                     'D=0|FA=3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,'
                                     '3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3'
                                     ',3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,'
                                     '3,3,3,3,3,3,3,3,3,3,3,3,3,3,3|PV=101,10'
                                     '2,102,101,102,114,111,114,111,113,106,1'
                                     '13,108,113,114,114,112,113,113,111,113,'
                                     '87,109,113,112,113,113,113,114,111,114,'
                                     '113,114,111,113,113,114,113,114,106,113'
                                     ',112,114,113,103,113,111,113,114,114,11'
                                     '4,112,110,111,113,113,114,114,110,112,1'
                                     '13,113,113,114,114,112,114,113,108,108'
                                     '\n',
                                     'AAATCGGGTCACTCCCACCTGAATACTGCGCTTTTCCGA'
                                     'TCGGCTTAAAAAATGGCGCACCACGAGATTA\n',
                                     '+\n', '~~~~~~~~~~~~~~~~~~~~~x~~~~~~~~'
                                     '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
                                     '~\n']
        conFq.close()
        self.filenames.append("TestR1.cons.fastq")

    def test_famsizestats(self):
        """
        TODO: Test GetFamSizeStats
        """
        pass

if __name__ == '__main__':
    unittest.main()
