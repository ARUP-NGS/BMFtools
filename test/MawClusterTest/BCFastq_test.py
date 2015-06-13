#!/user/bin/env python
import unittest
import sys
import subprocess
from utilBMF.HTSUtils import pFastqFile
from MawCluster.BCFastq import (pCompareFqRecsFast as cCompareFqRecsFast,
                                pairedFastqConsolidate, singleFastqConsolidate)

__author__ = 'dnephi'


class MyTestCase(unittest.TestCase):
    """
    TODO: Add a unit test to assert that the sum of the sizes of merged
    families equals the number of pre-merged reads.
    """

    def tearDown(self):
        pass

    def setUp(self):
        self.handle = pFastqFile("../data/TestR1.fastq")
        self.prefastqs = [i for i in self.handle]

    def test_dmp(self):
        assert cCompareFqRecsFast(self.prefastqs) == (
            '@MISEQ-M00736:68:000000000-A8D2D:1:2117:26553:9909 1:N:0:ACAGTG|'
            'FP=IndexPass|BS=AAAAAATGGACCCATTAACC|FM=3|ND=0|FA=3,3,3,3,3,3,3,'
            '3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,'
            '3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3|PV'
            '=101,102,102,101,102,114,111,114,111,113,106,113,108,113,114,114'
            ',112,113,113,111,113,87,109,113,112,113,113,113,114,111,114,113,'
            '114,111,113,113,114,113,114,106,113,112,114,113,103,113,111,113,'
            '114,114,114,112,110,111,113,113,114,114,110,112,113,113,113,114,'
            '114,112,114,113,108,108\nAAATCGGGTCACTCCCACCTGAATACTGCGCTTTTCCGA'
            'TCGGCTTAAAAAATGGCGCACCACGAGATTA\n+\n~~~~~~~~~~~~~~~~~~~~~x~~~~~~'
            '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')

    def test_pfc(self):
        pass

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
        subprocess.check_call("rm TestR1.cons.fastq", shell=True)
        conFq.close()

if __name__ == '__main__':
    unittest.main()
