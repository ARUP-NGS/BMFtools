#!/user/bin/env python
import shlex
import subprocess
import sys
import unittest

from BMFMain.Workflow import pairedFastqShades
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
        for filename in self.filenames:
            subprocess.check_call(["rm", filename])

    def setUp(self):
        self.filenames = []

    def test_compareCOBam(self):
        """
        TODO: Dummy case for making sure everything happens as it should.
        Check the following:
            1. Make sure that the is_qcfail flag is being set for reads
            with a run of 14 "A"s in a row in the barcode.
            2. Make sure that the is_qcfail flag is being set for reads
            with an N in the barcode.
            3. Make sure that the BS tag matches the read's name.
            4. Manually check that the PV tag matches the quality string.
            I would do
            assert sum(np.array(read.opt("PV").split(","), dtype=np.int32) == read.query_qualities) == len(read.query_qualities)
            
        """
        pass

if __name__ == '__main__':
    unittest.main()
