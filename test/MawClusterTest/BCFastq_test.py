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
    '''
    TODO: Add a unit test to assert that the sum of the sizes of merged
    families equals the number of pre-merged reads.
    '''

    def tearDown(self):
        for filename in self.filenames:
            subprocess.check_call(['rm', filename])

    def setUp(self):
        self.fq1 = '../data/lambdaTest_R1.fastq.gz'
        self.fq2 = '../data/lambdaTest_R3.fastq.gz'
        self.fqi = '../data/lambdaTest_R2.fastq.gz'
        self.filenames = []

    def test_FqProc(self):
        'What is this supposed to test?'
        PASS = True
        if not PASS:
            raise AssertionError(
                "Str: '%s' does not match unit test. Abort!" % tmpStr)

    def test_pfc(self):
        pairedFastqShades(self.fq1, self.fq2, indexFq=self.fqi)
        cStrToCheck = ("cat lambdaTest_R1.shaded.BS.cons.fastq | paste "
                       "- - - - | cut -d'|' -f4 | cut -d'=' -f2 | paste -sd+ "
                       "| bc")
        assert int(subprocess.check_output(cStrToCheck,
                                           shell=True).strip()) == 10000
        self.filenames += ['lambdaTest_R1.fastq.famstats.txt',
                           'lambdaTest_R1.shaded.BS.cons.fastq',
                           'lambdaTest_R1.shaded.BS.fastq',
                           'lambdaTest_R1.shaded.fastq',
                           'lambdaTest_R3.shaded.BS.cons.fastq',
                           'lambdaTest_R3.shaded.BS.fastq',
                           'lambdaTest_R3.shaded.fastq']

    def test_sfc(self):
        singleFastqConsolidate(self.fq1)
        with open('lambdaTest_R1.fastq.cons.fastq') as t:
            head = [next(t) for x in xrange(4)]
        assert head == ['@ 1:N:0:|FM=10000|ND=1082300|PV=3114,3114,3114,3114,3'
                        '114,3114,3114,3114,3114,3114,3114,3114,3114,3114,3114'
                        ',3114,3114,3114,3114,3114,3114,3114,3114,3114,3114,31'
                        '14,3114,3114,3114,3114,3114,3114,3114,3114,3114,3114,'
                        '3114,3114,3114,3114,3114,3114,3114,3114,3114,3114,311'
                        '4,3114,3114,3114,3114,3114,3114,3114,3114,3114,3114,3'
                        '114,3114,3114,3114,3114,3114,3114,3114,3114,3114,3114'
                        ',3114,3114,3114,3114,3114,3114,3114,3114,3114,3114,31'
                        '14,3114,3114,3114,3114,3114,3114,3114,3114,3114,3114,'
                        '3114,3114,3114,3114,3114,3114,3114,3114,3114,3114,311'
                        '4,3114,3114,3114,3114,3114,3114,3114,3114,3114,3114,3'
                        '114,3114,3114,3114,3114,3114,3114,3114,3114,3114,3114'
                        ',3114,3114,3114,3114,3114,3114,3114,3114,3114,3114,31'
                        '14,3114,3114,3114,3114,3114,3114,3114,3114,3114,3114,'
                        '3114,3114,3114,3114,3114|FA=2813,2821,2707,2590,2655,'
                        '2693,2795,2560,2555,2783,2567,2632,2632,2633,2611,262'
                        '4,2629,2657,2657,2734,2641,2695,2685,2571,2630,2583,2'
                        '623,2613,2598,2689,2598,2675,2664,2496,2577,2627,2569'
                        ',2562,2537,2586,2600,2547,2563,2580,2571,2549,2613,25'
                        '56,2700,2698,2638,2629,2613,2610,2601,2667,2765,2629,'
                        '2575,2617,2599,2563,2584,2622,2705,2637,2559,2713,259'
                        '2,2625,2671,2671,2615,2731,2774,2817,2781,2598,2702,2'
                        '746,2667,2728,2677,2622,2545,2641,2686,2624,2672,2604'
                        ',2594,2633,2602,2634,2649,2631,2597,2529,2586,2601,25'
                        '47,2587,2592,2552,2598,2545,2595,2596,2634,2576,2619,'
                        '2609,2712,2561,2677,2707,2644,2673,2652,2716,2677,261'
                        '2,2760,2616,2586,2705,2667,2738,2674,2710,2641,2654,2'
                        '644,2620,2579,2556,2737,2647,2581,2509,2667,2621,2548'
                        ',2648,2611,2751,2834\n', 'GATAGGATTATCACACGTCTGAACTCC'
                        'ATTCACTATTAAAACAAGATTAATCTCGTATAAAGTCTTCTGCTTAAAAAAAA'
                        'AAAATAAAAACAAAACATCTCACACCACACACATAAAAAAAAAAAAAAAAAAA'
                        'ACTAACTAATCAAA\n', '+\n', 'CBCCCCDDEFFCGGGGGGGGGGGHHH'
                        'HHHHHDHHHHHHHHHHHHHHHHGGHGHHHHHHHHHHHHHHHGHGHGGGGGGGG'
                        'GHHHHHGGEGGHHHHHHHHHGHGGGHHHHHHHHHHGHHGGGHGHFHHHHHHHG'
                        'HHHGHHHHHHHHHHH\n']
        self.filenames.append('lambdaTest_R1.fastq.cons.fastq')

    def test_famsizestats(self):
        '''
        TODO: Test GetFamSizeStats
        '''
        pass

if __name__ == '__main__':
    unittest.main()
