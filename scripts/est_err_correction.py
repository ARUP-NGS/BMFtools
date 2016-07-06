#!/usr/bin/env python
import sys
import pysam


def make_famstats_dict(path):
    ret = {}
    fh = open(path)
    line = next(fh)
    while("#Family size" not in line):
        line = next(fh)
    for line in fh:
        if line[0] == "#":
            break
        toks = line.strip().split()
        ret[int(toks[0])] = int(toks[1])
    return ret


def make_err_dict(path):
    return {int(line.strip().split()[0]): (float(line.strip().split()[1]),
                                           float(line.strip().split()[2])) for
            line in open(path) if line[0] != "#"}


def get_mean_err_correction(famstats_dict, err_dict, minFM=0):
    assert list(famstats_dict.keys()) == list(err_dict.keys())
    nfams = 0
    er_sum1, er_sum2 = 0., 0.
    for key in err_dict.keys():
        if key < minFM: continue
        fm = famstats_dict[key]
        nfams += fm
        err1, err2 = err_dict[key]
        er_sum1 += fm * err1
        er_sum2 += fm * err2
    return (er_sum1 / nfams, er_sum2 / nfams)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.stderr.write("Usage: %s errfm.txt famstats.txt <minFM>\n (minFM optional)\n" % sys.argv[0])
        sys.exit(1)
    maxFM = int(sys.argv[3]) if len(sys.argv) >= 4 else 0
    fm = make_famstats_dict(sys.argv[2])
    err = make_err_dict(sys.argv[1])
    sys.stdout.write("##%s Error rates.\n#minFM\tRead1\tRead2\n" % sys.argv[1])
    for m in range(maxFM):
        if m not in err: continue
        mean_err = get_mean_err_correction(fm, err, minFM=m)
        sys.stdout.write("%i\t%f\t%f\n" % tuple([m] + list(mean_err)))
    # Main goes here.
    sys.exit(0)
