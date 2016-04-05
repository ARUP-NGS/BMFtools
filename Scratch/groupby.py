#!/usr/bin/env python
import sys
import pysam
import itertools

def group_key(read):
    '''
    Note: the bam is already sorted by r1/r2. It's a little unsafe to do it this way,
    but it should still work.
    '''
    if read.is_read1:
        return (read.opt("SU"), read.opt("MU"),
                read.tid, read.next_reference_id,
                read.is_reverse, read.mate_is_reverse,
                read.is_read1, read.is_read1)
    else:
        return (read.opt("MU"), read.opt("SU"),
                read.next_reference_id, read.tid,
                read.mate_is_reverse, read.is_reverse,
                read.is_read2, read.is_read1)

def pysam_iterator(alignmentfile):
    for read in alignmentfile:
        if read.is_unmapped or read.mate_is_unmapped:
            continue
        yield read

def main():
    if len(sys.argv) < 2:
        print("Need a bam as positional argument 1.")
        sys.exit(1)
    r1dict = {}
    r2dict = {}
    for key, giter in itertools.groupby(pysam_iterator(pysam.AlignmentFile(sys.argv[1])), key=group_key):
        if key[7]:
            r1dict[":".join(map(str, key))] = list(giter)
        else:
            r2dict[":".join(map(str, key))] = list(giter)
    return {"r1": r1dict, "r2": r2dict}


if __name__ == "__main__":
    sys.exit(main())
