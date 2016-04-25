#!/usr/bin/env python
import pysam
import sys

if __name__ == "__main__":
    fmsum = 0
    rvsum = 0
    if len(sys.argv) < 2:
        sys.stderr.write("Usage: python %s <bampath> [2> results.txt]\n" % sys.argv[0].split("/")[-1])
        sys.exit(1)
    for read in pysam.AlignmentFile(sys.argv[1]):
        if read.flag & 2432:
            continue
        fmsum += read.opt("FM")
        rvsum += read.opt("RV")
    sys.stderr.write("%s-FM:%i;RV:%i.\n" % (sys.argv[1], fmsum, rvsum))
    sys.exit(0)
