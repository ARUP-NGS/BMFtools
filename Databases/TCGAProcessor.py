# TCGA Processor
# ####Components
# ####1. .maf --> .bed converter
# ####2. Mutation Frequency Calculator
# ####3. Mutation-Bedtool Mapper Across Samples
# ####
# ####Input:	InputGeneList (\n-delimited, HUGO naming)
# ####		EDIT: Not a .maf file
# ####		Set of .bed files with which to intersect the TCGA mutations
# ####
# ####Output:	(Tentative) List all mutations for each gene
# ####		List fraction of cases with a mutation in this gene
# ####
# ####		For each mutation type, list frequency (absolute number)
# ####		relative (fraction of samples with an alteration in this gene which have this specific alteration)
# ####
# ####		I'll figure out the rest of the output later
def main():
    import pybedtools
    import sys
    from pybedtools import BedTool
    import math
    import matplotlib
    matplotlib.use('Agg')
    import numpy as np
    from matplotlib import pyplot as plt
    from matplotlib_venn import venn3, venn3_unweighted, venn3_circles
    from stackedBarGraph import StackedBarGrapher
    	#Call Function to turn .maf to .bed {optional, not written yet 2/20/2014}
    SBG = StackedBarGrapher()

    chroms = [line.strip() for line in open(sys.argv[1] + '.1.temp')]
    startPos = [line.strip() for line in open(sys.argv[1] + '.2.temp')]
    endPos = [line.strip() for line in open(sys.argv[1] + '.3.temp')]
    geneName = [line.strip() for line in open(sys.argv[1] + '.4.temp')]
    vClass = [line.strip() for line in open(sys.argv[1] + '.5.temp')]
    vType = [line.strip() for line in open(sys.argv[1] + '.6.temp')]
    pID = [line.strip() for line in open(sys.argv[1] + '.7.temp')]
    FullLocs = []
    for i in range(len(chroms)):
    	FullLocs.append(str(chroms[i]) + ";" + str(startPos[i]) + ";" + str(endPos[i]) + ";" + geneName[i] + ";" + pID[i])
    uniquePinStudy = list(set(pID))
    totalPatients = len(uniquePinStudy)
    #Import Gene List
    geneList = [line.strip() for line in open(sys.argv[1]+'.gl')]
    dataBGlabels = []
    data_heights = []
    mappedPIDsets= []
    #Write header file for output
    bedRefList = [line.strip() for line in open(sys.argv[1] + ".bl")]
    outFile = open(str(sys.argv[1] + ".out.txt"), 'w')
    missedFile = open(str(sys.argv[1] + ".missed.txt"), 'w')
    missedFile.write("chr\tStart\tEnd\tHUGO_Alias\n")
    incCounts = np.zeros( (len(bedRefList) + 1, len(geneList)) )
    for l, beds in enumerate(bedRefList):
        coveredCounts = []
        geneCounts = []
        geneIndices = []
        bedRef = BedTool(beds)
        outFile.write("Covered\tMissed\tTotal\tGene Name\t {}\n".format(beds.split('/')[-1]))
        for gene in geneList:
            PtListForGene = []
            tempStr = ""
            for i, j in enumerate(geneName):
                if(j.lower() in gene.lower()): #changed == to 'in'
                    PtListForGene.append(pID[i])
                    tempStr += str(i)+";"
            geneIndices.append(tempStr)
            geneCounts.append(len(list(set(PtListForGene))))
	for i, gene in enumerate(geneList):
		fullInGene = []
		geneIndexTemp = geneIndices[i].split(";")[0:-1]
		for tmp in geneIndexTemp:
			fullInGene.append(FullLocs[int(tmp)])
		mappedPIDs = []
		for j, loc in enumerate(fullInGene):
			ersatzBF = BedTool("{} {} {} {}".format(loc.split(";")[0],loc.split(";")[1],loc.split(";")[2],loc.split(";")[3],loc.split(";")[4]),from_string=True)
			if(len(ersatzBF.intersect(bedRef)) > 0):
				mappedPIDs.append(loc.split(";")[4])
			elif l == 1:
				missedFile.write(str(ersatzBF))
		        ersatzBF.delete_temporary_history(ask=False)
                coveredCounts.append(len(list(set(mappedPIDs))))
		tmpP = ""
		for p in list(set(mappedPIDs)):
			tmpP += p + ";"
		mappedPIDsets.append(tmpP)
		#Time To Prepare Output - Human Readable
		outFile.write("{}\t{}\t{}\t{}\n".format(coveredCounts[-1], (geneCounts[i]-coveredCounts[-1]), geneCounts[i], geneList[i]))
	incCounts[l] = coveredCounts
    incCounts[-1] = geneCounts
    #Now I am trying to produce the stats for A-B, B-A, and A+B

        

    geneLabels = geneList
    for cnt, g in enumerate(geneLabels):
    	geneLabels[cnt] += " (Total Patients: " + str(geneCounts[cnt]) + ")"
    '''
    SBG.stackedBarPlot(ax, dataForBarGraph, d_colors, widths = d_widths, xLabels=geneLabels, yTicks=6, scale=True)
    plt.legend( ('Submitted','Designed','Total'))
    ax.set_xticklabels(geneLabels, rotation=18, multialignment = 'left', fontsize = 6 )
    plt.grid(True, which='minor')
    plt.title("Fraction of Patients in TCGA Covered for {}".format(str(sys.argv[1])), fontsize = 9)
    plt.ylabel("% Covered Patients With Mutation ", fontsize = 9)
    plt.tight_layout()
    plt.savefig(sys.argv[1]+'SumFreqNormalized')
    plt.close(fig)
    del fig/
    '''
    fig, ax = plt.subplots()
    width = 0.35
    ind = np.arange(len(geneCounts))
    plt.tight_layout()
    #Dummy Code
    plt.xticks(ind, geneLabels, rotation = 75, fontsize = 6)
    rects0 = ax.bar(ind, 100*incCounts[0]/geneCounts, width, color = "#FFFD00",edgecolor="k")
    rects1 = ax.bar(ind + width , 100*incCounts[1]/geneCounts, width, color = "#34D7A3",edgecolor="k")
    ax.set_title("Fraction of Patients in TCGA Covered for {}".format(str(sys.argv[1])), fontsize = 11)
    ax.set_ylabel("% Covered Patients With Mutation", fontsize = 8)
    ax.legend( ('Designed For Lung v1', 'Designed For Lung v2' ))
    ax.set_xlim(0,len(geneCounts)+width)
    def autolabel(rects):
    # Attaches text labels at location in barchart
	for cr, rect in enumerate(rects):
            if(not math.isnan(rect.get_x()) and not math.isnan(rect.get_height())):
                ax.text(rect.get_x() + rect.get_width()/2, 0.9*rect.get_height(), '%d'%int(rect.get_height() /100 * geneCounts[cr]), ha = 'center', va = 'bottom', fontsize = 8)
    autolabel(rects0)
    autolabel(rects1)
    #autolabel(rects2)
    #temp1 = []
    #for rect in rects0:
    #    temp1.append(rect.get_x())
    plt.tight_layout()
    plt.savefig(sys.argv[1] + 'SumCount.png')
    plt.close(fig)
    outFile.close()
    missedFile.close()

'''
	Tot = len(mutPtsUnion)
		outFile.close()
		#Print scaled image
		d_widths = [1.15] * len(geneList)
		#dataBGlabels contains the labels for these genes
		d_colors = ['#2166ac', '#fee090', '#fdbb84', '#fc8d59', '#e34a33', '#b30000']
		fig = plt.figure()
		ax =  fig.add_subplot(111)
		SBG.stackedBarPlot(ax, dataForBarGraph, d_colors, heights=data_heights, xLabels=dataBGlabels, yTicks=6, widths=d_widths, scale=False)
		plt.title("Summary for Frequencies of Given Mutations in TCGA")
		plt.tight_layout()
		plt.savefig(sys.argv[1]+'SumFreq')
		plt.close(fig)
		del fig
		
		fig = plt.figure()
		ax =  fig.add_subplot(111)
		SBG.stackedBarPlot(ax, dataForBarGraph, d_colors, xLabels=dataBGlabels, yTicks=6, widths=d_widths, scale=True)
		plt.title("Summary for Frequencies of Given Mutations in TCGA")
		plt.tight_layout()
		plt.savefig(sys.argv[1]+'SumFreqNormalized')
		plt.close(fig)
		del fig
'''
	#Program Complete
if __name__ == "__main__":
    main()

