#!/usr/bin/env python

def rev_test(inr, outr, minFrac=0.0, minPV=0, cap=93, minFM=0):
    if inr.opt("FM") < minFM:
        return 0
    # They have the same family size.
    for qin, pvin, fain, qout, pvout, faout in zip(inr.query_qualities, inr.opt("PV")[::-1], inr.opt("FA")[::-1],
            outr.query_qualities, outr.opt("PV")[::-1], outr.opt("FA")[::-1]):
        if pvin < minPV or fain * 1. / inr.opt("FM") < minFrac:
            if qout != 2:
                print("qout not 2 as expected (%i) in %s" % (qout, inr.query_name))
                return 0
        else:
            if qout != cap:
                print("qout not cap (%i) as expected (%i) in %s" % (cap, qout, inr.query_name))
                return 0
    return 1

def fw_test(inr, outr, minFrac=0.0, minPV=0, cap=93, minFM=0):
    if inr.opt("FM") < minFM:
        return 0
    # They have the same family size.
    for qin, pvin, fain, qout, pvout, faout in zip(inr.query_qualities, inr.opt("PV"), inr.opt("FA"),
            outr.query_qualities, outr.opt("PV"), outr.opt("FA")):
        if pvin < minPV or fain * 1. / inr.opt("FM") < minFrac:
            if qout != 2:
                print("qout not 2 as expected (%i) in %s" % (qout, inr.query_name))
                return 0
        else:
            if qout != cap:
                print("qout not cap (%i) as expected (%i) in %s" % (cap, qout, inr.query_name))
                return 0
    return 1



if __name__ == "__main__":

    import unittest
    import sys
    from subprocess import check_output
    import pysam
    import os.path

    # inbam = "FINAL_S376x3_Unpurified_GCATAACG-GCATAACG_L001_R1_001.R1.rsq.bam"
    inbam = "bigfams.bam"
    outbam = "capped.bam"
    minFrac = 0.7
    minFM = 2
    cap = 49
    minPV = 137
    check_output("bmftools cap -f%f -m%i -c%i -t%i %s %s" % (minFrac, minFM, minPV, cap, inbam, outbam), shell=True)
    class Test(unittest.TestCase):
        def test_cap(self):
            for inr, outr in zip(pysam.AlignmentFile(inbam), pysam.AlignmentFile(outbam)):
                if inr.is_reverse:
                    self.assertTrue(rev_test(inr, outr, minFrac=minFrac, minFM=minFM, cap=cap, minPV=minPV))
                else:
                    self.assertTrue(fw_test(inr, outr, minFrac=minFrac, minFM=minFM, cap=cap, minPV=minPV))


    unittest.main()
