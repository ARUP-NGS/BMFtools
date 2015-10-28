import numpy as np
import matplotlib.pyplot as plt
import argparse
import csv
import matplotlib
from matplotlib import gridspec
matplotlib.style.use('ggplot')
font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 22}

matplotlib.rc('font', **font)


def get_args():
    parser = argparse.ArgumentParser(description="Family Size Histogram")
    parser.add_argument("famstats", help="output from BMFTools famstats")
    parser.add_argument("binSize", help="size of bin for barplot", type=int)
    parser.add_argument("output", help="name of file to output to")
    return parser.parse_args()


def getBins(famstats, binSize):
    allValDict = {}
    lowValDict = {}
    currbin = 10
    with open(famstats) as t:
        for line in csv.reader(t,delimiter="\t"):
            if("#" in line[0]):
                continue
            size = int(line[0])
            count = int(line[1])
            if(size < 11):
                lowValDict[size] = count
            if(size < currbin+binSize+1):
                try:
                    allValDict[currbin] += count
                except KeyError:
                    allValDict[currbin] = count
            else:
                currbin += binSize
                allValDict[currbin] = count
    return allValDict, lowValDict

def famHistogram(famstats, binSize, output):
    highBinDict, lowBinDict = getBins(famstats, binSize)
    highFamSizes = np.array(np.zeros(len(highBinDict.keys())))
    highCounts = np.array(np.zeros(len(highBinDict.keys())))
    highRawCounts = np.array(np.zeros(len(highBinDict.keys())))
    lowFamSizes = np.array(np.zeros(len(lowBinDict.keys())))
    lowCounts = np.array(np.zeros(len(lowBinDict.keys())))
    lowRawCounts = np.array(np.zeros(len(lowBinDict.keys())))
    for index, key in enumerate(sorted(highBinDict.keys())):
        highFamSizes[index] = key
        highCounts[index] = np.log2(highBinDict[key])
        highRawCounts[index] = np.log2(highBinDict[key]*key)
    for index, key in enumerate(sorted(lowBinDict.keys())):
        lowFamSizes[index] = key
        lowCounts[index] = np.log2(lowBinDict[key])
        lowRawCounts[index] = np.log2(lowBinDict[key]*key)

    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 4])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1], sharey=ax1)
    ax3 = ax2.twinx()
    rects1 = ax1.bar(lowFamSizes, lowCounts, color='r')
    rects2 = ax2.bar(highFamSizes, highCounts, color='r',
                    width=max(highFamSizes)/(len(highFamSizes)+1))
    dense1 = ax3.plot(highFamSizes, highRawCounts, color='b', linewidth=5)
    ax1.set_xlim(1,11)
    ax2.set_xlim(11,max(highFamSizes))
    ax1.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax1.yaxis.tick_left()
    ax1.tick_params(labeltop='off')
    ax2.tick_params(left='off')
    ax2.tick_params(labelleft='off')
    ax1.axis('tight')
    ax2.axis('tight')
    ax2.set_xlabel("Family Size")
    ax1.set_ylabel("$Log_2$ observations of Family Sizes", color = 'r')
    ax3.set_ylabel("$Log_2$ total reads in families of that size", color = 'b')
    ax3.tick_params(axis='y', colors='blue')
    ax1.tick_params(axis='y', colors='red')
    plt.subplots_adjust(wspace=0.075)
    plt.savefig(output+".png")
    plt.show()



def main():
    args = get_args()
    famHistogram(args.famstats, args.binSize, args.output)
    return 0


if __name__ == "__main__":
    main()
