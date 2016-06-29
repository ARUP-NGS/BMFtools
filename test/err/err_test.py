#!/usr/bin/env python

import os
import subprocess
import sys
import hashlib

executables = ["bmftools", "bmftools_p", "bmftools_db"]

def main():
    try:
        genome_path = sys.argv[1]
    except IndexError:
        #raise ValueError("genome path (python err_test.py <genome_path>) required.")
        genome_path = "no_index_provided"
    if(os.path.isfile(genome_path) is False):
        print("WARNING: to test error rates, provide a path to build 37 human"
             " reference. This test passes automatically if not found.")
        return 0
    for ex in executables:
        subprocess.check_call("%s err fm -o err_test.out "
                              "%s NA12878.on_target.bam" % (ex, genome_path),
                              shell=True)
        if sys.version_info.major >= 3:
            assert "f4257518cfaccb82085f4d7708222803" == hashlib.md5(open("err_test.out", "r").read().encode()).hexdigest()
        else:
            assert "f4257518cfaccb82085f4d7708222803" == hashlib.md5(open("err_test.out", "r").read()).hexdigest()
    return 0

if __name__ == "__main__":
    sys.exit(main())
