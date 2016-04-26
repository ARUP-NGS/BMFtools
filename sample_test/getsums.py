#!/usr/bin/env python
import pysam
import sys

if __name__ == "__main__":
    fmsum = 0
    rvsum = 0
    ex = sys.argv[0].split("/")[-1]
    if len(sys.argv) < 2:
        sys.stderr.write("Usage: python %s <bampath> <outfile> "
                         "[omit outfile to write to stdout]\n" % ex)
        sys.exit(1)
    handle = open(sys.argv[2], "w") if len(sys.argv >= 3) else sys.stdout
    for read in pysam.AlignmentFile(sys.argv[1]):
        # 2432 = (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FREAD2)
        if read.flag & 2432:
            continue
        fmsum += read.opt("FM")
        rvsum += read.opt("RV")
    handle.write("%s-FM:%i;RV:%i.\n" % (sys.argv[1], fmsum, rvsum))
    sys.exit(0)
