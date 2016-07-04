#!/usr/bin/env python
import sys
import pysam
import operator
try:
    from cytoolz import frequencies as freq
except ImportError:
    from collections import Counter as freq

class FreqMaster:

    def set_position(self, contig, position):
        if not (contig and position):
            sys.stderr.write("No position provided. Doing nothing.\n")
            self.pci = self.pc = None
            return None
        if isinstance(contig, int):
            contig = self.fp.references[contig]
        self.pci = self.fp.pileup(contig, position)
        self.pc = next(self.pci)
        while self.pc.pos < position:
            self.pc = next(self.pci)

    def fill_tables(self):
        self.nucs = set(pr.alignment.seq[pr.query_position] for
                        pr in self.pc.pileups)
        self.tables = {nuc: freq(pr.alignment.opt("FM") for pr in
                                 self.pc.pileups if
                                 pr.alignment.seq[pr.query_position] == nuc)
                       for nuc in self.nucs}

    def __init__(self, bampath, contig=None, position=None):
        self.fp = pysam.AlignmentFile(bampath)
        self.tables = None
        self.set_position(contig, position)
        assert self.pc.pos == position
        self.fill_tables()

    def write_tables(self, outpath=None):
        if not self.tables:
            raise ValueError("Need tables to write them!")
        outhandle = open(outpath, "w") if outpath else sys.stdout
        outhandle.write("##Family Size\tFrequency\n")
        for nuc in sorted(self.nucs):
            outhandle.write("#Allele %c Frequency\n" % nuc)
            for fm, count in sorted(self.tables[nuc].items(),
                                    key=operator.itemgetter(0)):
                outhandle.write("%i\t%i\n" % (fm, count))
        outhandle.close()

if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.stderr.write("Usage: %s in.bam contig_name "
                         "zero_based_pos <outpath>\n"
                         "Omitting the outpath causes the program "
                         "to emit to stdout.\n" % sys.argv[0])
        sys.exit(1)
    table_writer = FreqMaster(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    table_writer.fill_tables()
    table_writer.write_tables(None if len(sys.argv) == 4 else sys.argv[4])
    sys.exit(0)
