import pysam
import sys

def freq(iterable):
    ret = {}
    for item in iterable:
        if item in ret:
            ret[item] += 1
        else:
            ret[item] = 1
    return ret


DEFAULT_PADDING = 1500
DEFAULT_INSERT_MAX = 500

def get_nearby(bampath, contig, start, stop, padding=DEFAULT_PADDING, insert_max=DEFAULT_INSERT_MAX):
    return [i for i in pysam.AlignmentFile(bampath).fetch(contig, start - padding, stop + padding)
            if (i.flag & 1796) == 0 and i.tlen <= insert_max]


class ReadPair(object):
    def __init__(self, r1, r2):
        self.r1 = r1
        self.r2 = r2

def segregate(bam, reclist):
    assert isinstance(bam, pysam.calignmentfile.AlignmentFile)
    names = freq(rec.qname for rec in reclist)
    singletons = []
    pair_hash = {}
    for rec in reclist:
        if names[rec.qname] == 1:
            singletons.append(rec)
        elif names[rec.qname] != 2:
            raise ValueError("This should be 2!")
        elif rec.qname in pair_hash:
            pair_hash[rec.qname] = ReadPair(pair_hash[rec.qname], rec)
        else:
            pair_hash[rec.qname] = rec
    return singletons, list(pair_hash.itervalues())

def get_olap(read, start, stop):
    if read.is_unmapped:
        return 0
    return sum(pos >= start and pos < stop for pos in read.get_reference_positions())

def pair_spans(pair, start, stop):
    return ((pair.r2.aend <= start and pair.r1.pos >= stop) or
            (pair.r1.aend <= start and pair.r2.pos >= stop))

def max_olap(pair, start, stop):
    return max((get_olap(pair.r1, start, stop), (get_olap(pair.r2, start, stop))))


def process(bampath, contig, start, stop, padding=DEFAULT_PADDING, insert_max=DEFAULT_INSERT_MAX,
            outpath=None, histpath=None):
    bam = pysam.AlignmentFile(bampath)
    singletons, pairs = segregate(bam, get_nearby(bampath, contig, start, stop, padding=padding, insert_max=insert_max))
    n_single = len(singletons)
    n_pairs = len(pairs)
    total = n_single + 2 * n_pairs
    singletons = [i for i in singletons if get_olap(i, start, stop) == 0]
    pairs = [i for i in pairs if not max_olap(i, start, stop)]
    max_olap_hist = freq(max_olap(pair, start, stop) for pair in pairs)
    if not histpath:
        histpath = bampath + ".hist.txt"
    print("Writing olap histogram for pairs not spanning to file '%s'" % histpath)
    with open(histpath, "w") as handle:
        handle.write("#Max probe overlap\tNumber of pairs\n")
        handle.write("\n".join("%i\t%i" % (olap, count) for olap, count in max_olap_hist.iteritems()))
    pairs = [pair for pair in pairs if not pair_spans(pair, start, stop) and pair.r1.is_proper_pair]
    for pair in pairs:
        assert not pair_spans(pair, start, stop)
        assert not max_olap(pair, start, stop)
    total_outside = len(singletons) + 2 * len(pairs)
    if not outpath:
        outpath = bampath + ".outof.%s.%i.%i.bam" % (contig, start, stop)
    with pysam.AlignmentFile(outpath, "wb", template=bam) as outbam:
        [outbam.write(single) for single in singletons]
        [(outbam.write(pair.r1), outbam.write(pair.r2)) for pair in pairs]
    print("Total in region: %i" % total)
    print("Total not pulled by probes in region: %i" % total_outside)
    print("Fraction in region unexpected: %f" % (float(total_outside) / total))
    print("Writing unexpected reads to %s" % outpath)

def main():
    return 0

if __name__ == "__main__":
    sys.exit(main())
