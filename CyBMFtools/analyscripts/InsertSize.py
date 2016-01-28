#!/usr/bin/env python

import array
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy.misc import factorial

'''
Sample data
# Insert sizes. Use `grep ^IS | cut -f 2-` to extract this part. The columns are: insert size, pairs total, inward oriented pairs, outward oriented pairs, other pairs
0	641	0	3	638
1	0	0	0	0
2	25	0	11	14
3	16	0	11	5
4	9	0	8	1
5	10	0	7	3
6	6	0	5	1
7	5	0	4	1
8	9	0	9	0
'''

def poisson(k, lamb):
    return (lamb**k / factorial(k) * np.exp(lamb))

class InsertSizeStats(object):

    alpha = 0.5

    def __init__(self, path, n=1500, bin_size=10):
        self.n = n
        self.n_bins = self.n / bin_size
        insert_sizes = array.array("I")
        counts = array.array("I")
        inward_counts = array.array("I")
        outward_counts = array.array("I")
        other_counts = array.array("I")
        for line in open(path, "r"):
            if line[:2] != "IS":
                continue
            toks = line.strip().split("\t")[1:]
            if toks[0][0] == "#":
                continue
            insert_sizes.append(int(toks[0]))
            counts.append(int(toks[1]))
            inward_counts.append(int(toks[2]))
            outward_counts.append(int(toks[3]))
            other_counts.append(int(toks[4]))
        self.insert_sizes = np.array(insert_sizes, dtype=np.uint32)
        self.counts = np.array(counts, dtype=np.uint32)
        self.inward_counts = np.array(inward_counts, dtype=np.uint32)
        self.outward_counts = np.array(outward_counts, dtype=np.uint32)
        self.other_counts = np.array(other_counts, dtype=np.uint32)

    def graph(self, outfname):
        with PdfPages(outfname) as pdf:
            data = self.counts[:self.n]
            n, bins, patches = plt.hist(data, self.n_bins, "r-", normed=1, facecolor="green", alpha=self.alpha)
            mu, std = norm.fit(data)
            y = mlab.normpdf(bins, mu, std)
            l = plt.plot(bins, y, 'r--', linewidth=1)
            plt.xlabel("Insert Size")
            plt.ylabel("Counts (pairs)")
            plt.grid(True)
            pdf.savefig()

    def __len__(self):
        return len(self.insert_sizes)
