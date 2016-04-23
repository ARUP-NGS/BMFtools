import pysam
try:
    from cytoolz import frequencies as freq
except ImportError:
    try:
        from stl.util import freq
    except ImportError:
        from collections import Counter as freq

def plp_pos(bam, contig, pos):
    """
    bampath, contig (string), position
    """
    return c


def aaf_pos(bam, contig, pos, ref='C', alt='T', minFM=1):
    a = pysam.AlignmentFile(bam)
    b = a.pileup(contig, pos)
    c = b.next()
    while c.pos < pos:
        c = b.next()
    counts = freq(i.alignment.seq[i.query_position] for i in c.pileups
                  if i.alignment.opt("FM") >= minFM)
    try:
        return float(counts[alt]) / counts[ref]
    except KeyError:
        if ref not in counts:
            return -137.
        return 0.
