from collections import defaultdict as dd
import pysam


def get_locus_counts(path, contig="7", pos=55249070, ref="C", alt="T"):
    refc, altc = 0, 0
    af = pysam.AlignmentFile(path)
    b = af.pileup(contig, pos, max_depth=500000)
    c = next(b)
    while c.pos < pos: c = next(b)
    for read in c.pileups:
        if not read.query_position: continue
        if read.alignment.seq[read.query_position] == ref: refc += 1
        if read.alignment.seq[read.query_position] == alt: altc += 1
    return refc, altc

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        sys.stderr.write("Usage: python %s <bedpath>.\n"
                         "Emits mean doc and egfr to stdout.\n" % sys.argv[0])
    sys.stdout.write("Name\tRef\tAlt\t%Alt\n")
    for path in sys.argv[1:]:
        ref, alt = get_locus_counts(path)
        sys.stdout.write("%s\t%i\t%i\t%f\n" % (path, ref, alt,
                                               alt * 100. / (ref + alt)))
