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


def main():
    for ex in ["bmftools_db", "bmftools"]:
        cstr = "%s hashdmp -o hashdmp_test.out hashdmp_test.fq" % ex
        print "Command string: %s" % cstr
        subprocess.check_call(shlex.split(cstr))
        fqh  = pysam.FastqFile("hashdmp_test.out")
        r1 = fqh.next()
        tags1 = get_tags(r1)
        assert tags1["FM"] == 7
        assert round(tags1["NF"], 2) == 0.14
        assert tags1["RV"] == 2
        assert tags1["DR"]
        assert len(r1.name) == 16
        r1 = fqh.next()
        tags1 = get_tags(r1)
        assert tags1["FM"] == 1
        try:
            assert tags1["FP"] == 0
        except AssertionError:
            print(str(r1))
            raise AssertionError



    return 

if __name__ == "__main__":
    sys.exit(main())
