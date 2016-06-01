#!/usr/bin/env python

import pysam


def generate_passed_lines(vcfpath):
    for var in pysam.VariantFile(vcfpath):
        if "BMF_VET" in var.info:
            bp = var.info["BMF_VET"]
            if sum(bp) != bp[0]:
                yield var


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        sys.stderr.write("Usage: %s <vcfpath>.\n" % sys.argv[0])
    count = 0
    for line in generate_passed_lines(sys.argv[1]):
        count += 1
        sys.stdout.write(str(line))
    sys.stderr.write("Wrote %i passing lines.\n" % count)
