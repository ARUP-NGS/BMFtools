#!/usr/bin/env python
import shlex
import subprocess
import sys
import unittest

from utilBMF.QC import GetAllQCMetrics

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

    def test_GetAllQCMetrics(self):
        """
        Test GetAllQCMetrics
        """
        pass

if __name__ == '__main__':
    unittest.main()
