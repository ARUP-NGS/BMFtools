__author__ = 'dnephi'

import unittest
import sys
from utilBMF.HTSUtils import pFastqFile
from MawCluster.BCFastq import (cFRF_helper as compareFqRecsFast,
                                pairedFastqConsolidate, singleFastqConsolidate)

class MyTestCase(unittest.TestCase):

    def tearDown(self):
        pass

    def setUp(self):
        self.handle = pFastqFile("TestR1.fastq")
        self.prefastqs = [i for i in self.handle]
        #self.handle.refresh()
        #create dummy ucsc references try
        #full_gene
        #THIS WILL SHOW ONLY THE MERGING OF REGION 1

    def test_dmp(self):
        assert compareFqRecsFast(self.prefastqs) == '@MISEQ-M00736:68:000000000-A8D2D:1:2117:26553:9909 1:N:0:ACAGTG|FP=IndexPass|BS=AAAAAATGGACCCATTAACC|FM=3|ND=0|FA=3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3|PV=101,102,102,101,102,114,111,114,111,113,106,113,108,113,114,114,112,113,113,111,113,87,109,113,112,113,113,113,114,111,114,113,114,111,113,113,114,113,114,106,113,112,114,113,103,113,111,113,114,114,114,112,110,111,113,113,114,114,110,112,113,113,113,114,114,112,114,113,108,108\nAAATCGGGTCACTCCCACCTGAATACTGCGCTTTTCCGATCGGCTTAAAAAATGGCGCACCACGAGATTA\n+\n~~~~~~~~~~~~~~~~~~~~~x~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n'

    def test_pfc(self):
        pass

    def test_sfc(self):
        singleFastqConsolidate("TestR1.fastq")
        conFq = open("TestR1.cons.fastq").readlines()
        assert conFq == ['@AAAAAATGGACCCATTAACC 1:N:0:ACAGTG|FP=IndexPass|BS=AAAAAATGGACCCATTAACC|FM=3|ND=0|FA=3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3|PV=101,102,102,101,102,114,111,114,111,113,106,113,108,113,114,114,112,113,113,111,113,87,109,113,112,113,113,113,114,111,114,113,114,111,113,113,114,113,114,106,113,112,114,113,103,113,111,113,114,114,114,112,110,111,113,113,114,114,110,112,113,113,113,114,114,112,114,113,108,108\n', 'AAATCGGGTCACTCCCACCTGAATACTGCGCTTTTCCGATCGGCTTAAAAAATGGCGCACCACGAGATTA\n', '+\n', '~~~~~~~~~~~~~~~~~~~~~x~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n']

if __name__ == '__main__':
    unittest.main()
