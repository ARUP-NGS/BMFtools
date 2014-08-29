from Bio import SeqIO
import argparse
import BarcodeFastqTools

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--fq', help="Provide your fastq file(s). Note: if '--BAM' option is set, this needs to be the Family Fastq files. Note: must the the output of ", nargs = "+", metavar = ('reads'))
    parser.add_argument('-o','--outfq', help="Sets a destination file for the fastq. If not set, follows a default", metavar = ('reads'),default="default")
    parser.add_argument('-s','--stringency',help="Provide a stringency - fraction of family members who must entirely agree to be written.",default="0.9")
    args=parser.parse_args()
    outFq = args.outfq
    stringency = args.stringency
    inFq = SeqIO.parse(args.fq[0],'fastq')
    if(outFq=="default"):
        outFq= '.'.join(args.fq[0].split('.')[0:-1]) + 'cons.fastq'
    outputHandle = open(outFq,'w')
    workingBarcode = ""
    workingSet = []
    for fqRec in inFq:
        barcode4fq = fqRec.description.split("###")[-2].strip()
        #print("Working barcode: {}. Current barcode: {}.".format(workingBarcode,barcodeRecord))
        #print("name of read with this barcode: {}".format(record.qname))
        #print("Working set: {}".format(workingSet))
        if("AAAAAAAAAA" in barcode4fq or "TTTTTTTTTT" in barcode4fq or "CCCCCCCCCC" in barcode4fq or "GGGGGGGGGG" in barcode4fq):
            continue
        if(int(fqRec.description.split('###')[-1].strip()) < 2):
            continue
        if(workingBarcode == ""):
            workingBarcode = barcode4fq
            workingSet = []
            workingSet.append(fqRec)
            continue
        elif(workingBarcode == barcode4fq):
            workingSet.append(fqRec)
            continue
        elif(workingBarcode != barcode4fq):
            mergedRecord, success = BarcodeFastqTools.compareFastqRecords(workingSet,stringency=float(stringency))
            if(success == True):
                SeqIO.write(mergedRecord, outputHandle, "fastq")
            workingSet = []
            workingSet.append(fqRec)
            workingBarcode = barcode4fq
            continue
        else:
            raise RuntimeError("No idea what's going on. This code should be unreachable")
    inFq.close()
    outputHandle.close()
    return "SMILES"

if(__name__=="__main__"):
    main()