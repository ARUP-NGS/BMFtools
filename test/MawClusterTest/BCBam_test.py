#!/user/bin/env python
import shlex
import subprocess
import sys
import unittest
import pysam
import numpy as np
from BMFMain.Workflow import pairedFastqShades
from utilBMF.HTSUtils import pFastqFile, PipeAlignTag
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
        self.fq1 = "../data/lambdaTest_R1.shaded.BS.cons.fastq.gz"
        self.fq2 = "../data/lambdaTest_R3.shaded.BS.cons.fastq.gz"
        self.ref = "../data/OurPhageRef/OurPhageLambda.fasta"
        self.badBC = ['A'*16, 'G'*16, 'C'*16, 'T'*16]
        self.filenames = []

    def test_compareCOBam(self):
        PipeAlignTag(self.fq1, self.fq2, ref=self.ref)
        bam = pysam.AlignmentFile(
            "../data/lambdaTest_R1.shaded.BS.cons.fastq.mem.bam")
        self.filenames += ([
            "../data/lambdaTest_R1.shaded.BS.cons.fastq.mem.bam",
            "../data/lambdaTest_R1.shaded.BS.cons.fastq.mem.bam.bai",])
        for read in bam:
            if read.qname in self.badBC:
                if not read.is_qcfail:
                    raise AssertionError(
                            "BS: %s does not properly fail QC" % (read.qname))
            if 'N' in read.qname:
                if not read.is_qcfail:
                    raise AssertionError(
                            "BS: %s does not properly fail QC" % (read.qname))
            if read.qname != read.opt('BS'):
                raise AssertionError("Read: %s has mismatch in name and BS"
                                     " tag" % (read.qname))
            if read.opt('PV') != read.query_qualities:
                pv = [i if(i < 94) else(93) for i in  read.opt('PV')]
                if pv != list(read.query_qualities):
                    raise AssertionError("Read: %s has mismatched quality"
                                         " array and PV tag\nPV Tag:\n%s\n"
                                         "query_qualities: %s" % (
                                            read.qname,
                                            pv,
                                            read.query_qualities))

if __name__ == '__main__':
    unittest.main()
