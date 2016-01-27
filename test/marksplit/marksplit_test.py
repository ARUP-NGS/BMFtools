import sys
import subprocess
import shlex

import pysam
import numpy as np

def get_tags(read):
    tag_strs = read.comment.split("\t")[1:]
    ret = {}
    for tag in tag_strs:
        if tag[3] == 'i':
            ret[tag[:2]] = int(tag.split(":")[2])
        elif tag[3] == 'f':
            ret[tag[:2]] = float(tag.split(":")[2])
        elif tag[3] == 'B':
            ret[tag[:2]] = np.array(tag.split(":")[2].split(",")[1:], dtype=np.uint32)
        else:
            raise NotImplementedError("Haven't finshed this.")
    return ret

def check_bc(read):
    bc = read.comment.split("=")[2]
    fp = int(read.comment.split("|")[1].split("=")[1])
    assert bc[0] in "ZFR"
    if "N" in bc:
        assert fp == 0

def marksplit_main():
    for ex in ["bmftools_db", "bmftools", "bmftools_p"]:
        cstr = ("../../%s dmp -n0 -sTGACT -t12 -o marksplit_test_tmp -l 10 "
                "-v 11 marksplit_test.R1.fq marksplit_test.R2.fq" % ex)
        print "Command string: %s" % cstr
        subprocess.check_call(shlex.split(cstr))
        for read in pysam.FastqFile("marksplit_test_tmp.tmp.0.R1.fastq"):
            check_bc(read)
    return 0

def main():
    marksplit_main()

if __name__ == "__main__":
    sys.exit(main())
