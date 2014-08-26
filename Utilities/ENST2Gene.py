#/mounts/anaconda/bin/python

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-c','--counts', help="Provide your tab-delimited, three-column counts file.", metavar = ('counts'),required=True)
    parser.add_argument('-t','--table',help="Table containing ENST entries and corresponding gene names",required=True)
    args=parser.parse_args() 
    countsFile = args.counts
    table = args.table
    outputFile,countReplaced = ENST2Gene(countsFile,table)
    print("Output file is: {}".format(outputFile))
    print("Number of IDs replaced is {}".format(countReplaced))
    return

def ENST2Gene(countsFile,table,outfile="default"):
    import csv
    #Parse in the table
    tableFile = open(table,"r")
    dictEntries = [line.split() for line in tableFile]
    BarDict = {}
    for entry in dictEntries:
        BarDict[entry[0]]=entry[1]
    tableFile.close()
    counts = open(countsFile,"r")
    if(outfile=="default"):
        outfile='.'.join(countsFile.split('.')[0:-1])+'.mergedConv.counts.txt'
    write=open(outfile,"w",0)
    countsEntries = [line.split() for line in counts]
    SpamAndEggs = csv.writer(write,delimiter="\t")
    counter = 0
    for countRecord in countsEntries:
        try:
            newID = BarDict[countRecord[1]]
            newRecord = [countRecord[0],newID,countRecord[2]]
            counter+=1
        except KeyError:
            newRecord=countRecord
        except:
            newRecord=countRecord
        SpamAndEggs.writerow(newRecord)
    write.close()
    return outfile, counter

if(__name__=="__main__"):
    main()
