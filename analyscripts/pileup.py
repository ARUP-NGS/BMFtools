'''
Unexpected procotols...

To get the number of non-reference mapped bases from a VCF,
do these two things! :)

```bash
grep -v ^# $Infile| cut -f8 | grep -v 'INDEL' | \
cut -d';' -f2 | cut -d'=' -f2 | cut -d',' -f3,4 | sort | uniq -c
```
sum(int(k[0]) * sum(map(int, k[1].split(","))) for k in [
    i.strip().split() for i in
    open("NoBMF.nonref.hist", "r").read().split("\n") if i != '']
    if sum(map(int, k[1].split(","))) < 100)

If doing it with python after tabixing:
np.sum(np.array(list(cfi([
    dict(tup.split("=") for tup in
    i.info.split(";"))['I16'].split(",")[2:4] for i in
    inHandle if "INDEL" not in i.info])), dtype=np.int))

'''


def GetSumOfDifferencesFromTheReference(vcfpath):
    from subprocess import check_call
    from utilBMF.HTSUtils import TrimExt
    import pysam
    import numpy as np
    from sys import stderr
    from itertools import chain
    cfi = chain.from_iterable
    bgvcfpath = TrimExt(vcfpath) + ".gz"
    check_call("bgzip -c %s > %s" % (vcfpath, bgvcfpath), shell=True)
    stderr.write("bgvcf now at %s" % bgvcfpath)
    tabixstr = "tabix " + bgvcfpath
    stderr.write("Now calling tabixstr: '%s'" % tabixstr)
    check_call("tabix %s" % bgvcfpath, shell=True)
    infh = open(bgvcfpath, "rb")
    tabixhandle = pysam.tabix_iterator(infh, pysam.asVCF())
    return np.sum(np.array(list(cfi([dict(tup.split("=") for
                                          tup in i.info.split(";"))[
        'I16'].split(",")[2:4] for i in tabixhandle if
                                     "INDEL" not in i.info])), dtype=np.int64))


if __name__ == "__main__":
    import sys
    msg = ("Total non-reference aligned bases:"
           " %s\n" % GetSumOfDifferencesFromTheReference(sys.argv[-1]))
    sys.stdout.write(msg)
    sys.exit(0)
