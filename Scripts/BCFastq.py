from Bio import SeqIO
import pysam
from HTSUtils import printlog as pl


def BarcodeSort(inFastq, outFastq="default"):
    pl("Sorting {} by barcode sequence.".format(inFastq))
    from subprocess import call
    if(outFastq == "default"):
        outFastq = '.'.join(inFastq.split('.')[0:-1]) + '.BS.fastq'
    BS1 = "cat " + inFastq + " | paste - - - - | grep 'HomingPass' | sed "
    BS1 += "'s:###:###\t:g' |  awk 'BEGIN {{FS=OFS=\"\t\"}};{{print $3,$0}}' |"
    BS1 += " sort -k1,1 | awk 'BEGIN {{OFS=FS=\"\t\"}};"
    BS1 += "{{print $2,$3,$4,$5,$6,$7,$8,$9,$10}}' | "
    BS1 += "sed 's:###\t:###:g' | sed 's:\t$::g' "
    BS1 += "| sed 's:\t$::g' | tr '\t' '\n' > "
    BSstring = BS1 + outFastq
    call(BSstring, shell=True)
    pl("Command: {}".format(BSstring.replace(
        "\t", "\\t").replace("\n", "\\n")))
    return outFastq


def compareFastqRecords(R, stringency=0.9):
    from Bio.SeqRecord import SeqRecord
    seqs = [str(record.seq) for record in R]
    maxScore = 0
    Success = False
    for seq in seqs:
        # print("Seq: {}".format(str(seq)))
        numEq = sum(str(seq) == str(seqItem) for seqItem in seqs)
        if(numEq > maxScore):
            maxScore = numEq
            finalSeq = str(seq)
    frac = numEq * 1.0 / len(R)
    PASS = frac > stringency
    # print("Fraction {}. Stringency: {}. Pass? {}.".format(
    # frac,stringency,PASS))
    if(PASS):
        Success = True
    consolidatedRecord = SeqRecord(seq=finalSeq, id=R[0].id,
                                   letter_annotations=R[0].letter_annotations,
                                   name=R[0].name,
                                   description=R[0].description)
    return consolidatedRecord, Success


def HomingSeqLoc(fq, homing, bar_len=12):
    pl("Now beginning HomingSeqLoc.")
    InFastq = SeqIO.parse(fq, "fastq")
    Tpref = '.'.join(fq.split('.')[0:-1])
    Prefix = Tpref.split('/')[-1]
    StdFilename = Prefix + '.{}.fastq'.format("homing" + str(bar_len))
    ElseFilename = Prefix + '.else.fastq'
    ElseLoc = Prefix + '.else.supp'
    StdFastq = open(StdFilename, 'w', 0)  # Homing at expected Location
    ElseFastq = open(ElseFilename, 'w', 0)
    ElseLocations = open(ElseLoc, 'w', 0)
    for read in InFastq:
        seq = str(read.seq)
        if(seq.find(homing) == -1):
            read.description += " ###HomingFail"
            SeqIO.write(read, StdFastq, "fastq")
        elif(seq[bar_len:bar_len + len(homing)] == homing):
            read.description += " ###HomingPass"
            SeqIO.write(read, StdFastq, "fastq")
        else:
            read.description = " ###HomingFail"
            SeqIO.write(read, StdFastq, "fastq")
            ElseLocations.write(repr(seq.find(homing)) + "\t" +
                                read.name + "\n")
    StdFastq.close()
    ElseFastq.close()
    ElseLocations.close()
    return StdFilename, ElseFilename


def fastq_sort(in_fastq, out_fastq):
    pl("Now beginning fastq_sort.")
    import subprocess
    outfile = open(out_fastq, 'w')
    command_str = 'cat {} | paste - - - - | '.format(in_fastq)
    command_str += 'sort -k1,1 -t " " | tr "\t" "\n"'
    subprocess.call(command_str, stdout=outfile, shell=True)
    outfile.close()
    return(command_str)


def FastqRegex(fq, string, matchFile="default", missFile="default"):
    from subprocess import call
    if(matchFile == "default"):
        matchFile = ('.'.join(fq.split(
                     '.')[0:-1]) + '.match.fastq').split('/')[-1]
    if(missFile == "default"):
        missFile = (
            '.'.join(fq.split('.')[0:-1]) + '.miss.fastq').split('/')[-1]
    CommandStr = "cat {} | paste - - - - | grep ".format(fq)
    CommandStr += "'{}' | tr '\t' '\n' > {}".format(string, matchFile)
    call(CommandStr, shell=True)
    CommandStr2 = "cat {} | paste - - - - | grep ".format(fq)
    CommandStr2 += "-v '{}' | tr '\t' '\n' > {}".format(
        string,
        missFile)
    call(CommandStr2, shell=True)
    return(CommandStr, CommandStr2, matchFile, missFile)


def fastx_trim(infq, outfq, n):
    pl("Now beginning fastx_trimmer.")
    import subprocess
    command_str = ['fastx_trimmer', '-l', str(n), '-i', infq, '-o', outfq]
    pl(command_str)
    subprocess.call(command_str)
    return(command_str)


def findProperPairs(infq1, infq2, index1="default", index2="default",
                    outfq1="default", outfq2="default", outfqSingle="default"):
    """Assuming that the information here is sorted by barcode..."""
    import collections
    pl("Now beginning findProperPairs.")
    from Bio import SeqIO
    if(index1 == "default"):
        raise ValueError(
            "An index for the barcodes in reads1 set must be present.")
    if(index2 == "default"):
        raise ValueError(
            "An index for the barcodes in reads2 set must be present.")
    if(outfq1 == "default"):
        outfq1 = '.'.join(infq1.split('.')[0:-1]) + '.proper.fastq'
    if(outfq2 == "default"):
        outfq2 = '.'.join(infq2.split('.')[0:-1]) + '.proper.fastq'
    if(outfqSingle == "default"):
        outfqSingle = '.'.join(infq1.split('.')[0:-1]) + '.solo.fastq'
    prefix = '.'.join(infq1.split('.')[0])
    outfq1Handle = open(outfq1, 'w')
    outfq2Handle = open(outfq2, 'w')
    outfqSingleHandle = open(outfqSingle, 'w')
    indexHandle1 = open(index1, "r")
    indexHandle2 = open(index2, "r")
    barcodeList1 = [l.strip().split()[-1] for l in indexHandle1.readlines()]
    barcodeList2 = [l.strip().split()[-1] for l in indexHandle1.readlines()]
    indexHandle1.close()
    fullList = barcodeList1 + barcodeList2
    del barcodeList1
    del barcodeList2
    shared = [x for x, y in collections.Counter(fullList).items() if y > 1]
    infq1Handle = SeqIO.parse(infq1, "fastq")
    infq2Handle = SeqIO.parse(infq2, "fastq")
    BarDict = {}
    for entry in shared:
        BarDict[entry] = "SpanishInquisition"
    indexHandle2.seek(0)
    readNum = 0
    for read in infq2Handle:
        queryBC = read.description.split('###')[-2]
        try:
            temp = BarDict[queryBC]
            # This is simply checking whether or not
            # it has a mate with the same barcode sequence.
            read.name = "{}:ProperPairRd#{}".format(prefix, readNum)
            readNum += 1
            SeqIO.write(read, outfq2Handle, "fastq")
        except KeyError:
            SeqIO.write(read, outfqSingleHandle, "fastq")

    readNum = 0
    for read in infq1Handle:
        queryBC = read.description.split('###')[-2]
        try:
            temp = BarDict[queryBC]
            read.name = "{}:ProperPairRd#{}".format(prefix, readNum)
            readNum += 1
            SeqIO.write(read, outfq1Handle, "fastq")
        except KeyError:
            SeqIO.write(read, outfqSingleHandle, "fastq")
    outfq1Handle.close()
    outfq2Handle.close()
    infq1Handle.close()
    infq2Handle.close()
    indexHandle2.close()
    return outfq1, outfq2, outfqSingleHandle


def PairFastqBarcodeIndex(tags_file1, tags_file2, index_file="default"):
    pl("Now beginning GenerateFullFastqBarcodeIndex for {} and {}.".format(
        tags_file1, tags_file2))
    from subprocess import call
    if(index_file == "default"):
        index_file = '.'.join(tags_file1.split('.')[0:-1]) + ".barIdx"
    cmd = "cat {} {} | sed 's:###::g' | paste - - - - | awk ".format(
        tags_file1, tags_file2)
    cmd += "'BEGIN {{FS=\"\t\"}};{{print $2}}' | sort | uniq -c | awk 'BEGIN "
    cmd += "{{OFS=\"\t\"}};{{print $1,$2}}' | sort -k1,1n > {}".format(
        index_file)
    pl("CommandStr = {}".format(cmd.replace("\t", "\\t")))
    call(cmd, shell=True)
    return index_file


def GenerateOnePairFastqBarcodeIndex(tags_file, index_file="default"):
    pl("Now beginning GenerateFullFastqBarcodeIndex for {}.".format(tags_file))
    from subprocess import call
    if(index_file == "default"):
        index_file = '.'.join(tags_file.split('.')[0:-1]) + ".barIdx"
    cmd = "cat {} | sed 's:###::g' | paste - - - - | awk ".format(tags_file)
    cmd += "'{{print $4}}' | sort | uniq -c | awk 'BEGIN "
    cmd += "{{OFS=\"\t\"}};{{print $1,$2}}' > {}".format(index_file)
    pl("CommandStr = {}".format(cmd.replace("\t", "\\t")))
    call(cmd, shell=True)
    return index_file


def GenerateSingleBarcodeIndex(tags_file, index_file="default"):
    pl("Now beginning GenerateSingleBarcodeIndex for {}.".format(tags_file))
    from subprocess import call
    if(index_file == "default"):
        index_file = '.'.join(tags_file.split('.')[0:-1]) + ".barIdx"
    cmd = "cat {} | sed 's:###: ###:g' | ".format(tags_file)
    cmd += "paste - - - - | grep -v \"HomingFail\" | "
    cmd += "awk 'BEGIN {{FS=\"\t\";OFS=\"\t\"}};"
    cmd += "{{print $2}}' | sort | uniq -c | awk 'BEGIN {{OFS=\"\t\"}}"
    cmd += ";{{print $1,$2}}' > {}".format(index_file)
    pl("CommandStr = {}".format(cmd.replace("\t", "\\t")))
    call(cmd, shell=True)
    return index_file


def GetFamilySizeSingle(
        trimfq,
        BarcodeIndex,
        outfq="default",
        singlefq="default"):
    pl("Running GetFamilySizeSingle for {}.".format(trimfq))
    infq = SeqIO.parse(trimfq, "fastq")
    if(outfq == "default"):
        outfq = '.'.join(trimfq.split('.')[0:-1]) + ".fam.fastq"
    outfqBuffers = open(outfq, "w", 0)
    index = open(BarcodeIndex, "r")
    if(singlefq == "default"):
        singlefq = '.'.join(
            trimfq.split('.')[
                0:-1]) + ".lonely.hearts.club.band.fastq"
    singlefqBuffer = open(singlefq, "w", 0)
    TotalReads, ReadsWithFamilies = 0, 0
    dictEntries = [line.split() for line in index]
    BarDict = {}
    for entry in dictEntries:
        BarDict[entry[1]] = entry[0]
    for read in infq:
        # index.seek(0)
        TotalReads += 1
        # print("description is {}".format(read.description))
        # print("-1 index of description is {}".format(
        # read.description.split("###")[-1].strip()))
        # print("-2 index of description is {}".format(
        # read.description.split("###")[-2].strip()))
        readTag = read.description.split("###")[-1].strip()
        newRead = read
        # print("readTag is _" + readTag + "_")
        try:
            famSize = BarDict[readTag]
        except KeyError:
            famSize = 0
        newRead.description = read.description + " ###" + str(famSize)
        # print("famSize = _{}_".format(str(famSize)))
        # print("The value of this comparison to 1 is {}".format(
        # str(famSize=="1")))
        if(str(famSize) == "0"):
            continue
        if(str(famSize) == "1"):
            print("Hey, I found a singleton")
            SeqIO.write(newRead, singlefqBuffer, "fastq")
        else:
            print("Hey, I found a read with a family!")
            ReadsWithFamilies += 1
            SeqIO.write(newRead, outfqBuffers, "fastq")
    outfqBuffers.close()
    singlefqBuffer.close()
    return outfq, TotalReads, ReadsWithFamilies


def mergeBarcodes(fq1, fq2, out1="default", out2="default"):
    if(out1 == "default"):
        out1 = '.'.join(fq1.split('.')[0:-1]) + '.mergeBC.fastq'
    if(out2 == "default"):
        out2 = '.'.join(fq2.split('.')[0:-1]) + '.mergeBC.fastq'
    reads1 = SeqIO.parse(fq1, "fastq")
    reads2 = SeqIO.parse(fq2, "fastq")
    outFq1 = open(out1, "w")
    outFq2 = open(out2, "w")
    while True:
        try:
            read1 = reads1.next()
            read2 = reads2.next()
            read1Desc = read1.description.split('###')
            read2Desc = read2.description.split('###')
            if(read1Desc[1].strip() == "HomingFail"):
                continue
            if(read2Desc[1].strip() == "HomingFail"):
                continue
            newBarcode = read1Desc[-1] + read2Desc[-1]
            read1Desc[-1] = newBarcode
            read2Desc[-1] = newBarcode
            read1.description = "###".join(read1Desc)
            read2.description = "###".join(read2Desc)
            SeqIO.write(read1, outFq1, "fastq")
            SeqIO.write(read2, outFq2, "fastq")
        except StopIteration:
            break
    reads1.close()
    reads2.close()
    return out1, out2


def mergeSequencesFastq(fq1, fq2, output="default"):
    pl("mergeSequencesFastq for {} and {}".format(fq1, fq2))
    if(output == "default"):
        output = fq1.split('.')[0] + '.merged.fastq'
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    reads1 = SeqIO.parse(fq1, "fastq")
    reads2 = SeqIO.parse(fq2, "fastq")
    outFastq = open(output, "w", 0)
    while True:
        try:
            read1 = reads1.next()
            read2 = reads2.next()
            outread = read1
            outread = SeqRecord(
                Seq(str(read1.seq) + str(read2.seq), "fastq"),
                id=read1.id, description=read1.description)
            outread.letter_annotations[
                'phred_quality'] = read1.letter_annotations[
                'phred_quality'] + read2.letter_annotations['phred_quality']
            SeqIO.write(outread, outFastq, "fastq")
        except StopIteration:
            break
    reads1.close()
    reads2.close()
    outFastq.close()
    return output


def pairedFastqConsolidate(fq1, fq2, outFqPair1="default",
                           outFqPair2="default", stringency=0.9):
    pl("Now running pairedFastqConsolidate on {} and {}.".format(fq1, fq2))
    if(outFqPair1 == "default"):
        outFqPair1 = '.'.join(fq1.split('.')[0:-1]) + 'cons.fastq'
    if(outFqPair2 == "default"):
        outFqPair2 = '.'.join(fq2.split('.')[0:-1]) + 'cons.fastq'
    inFq1 = SeqIO.parse(fq1, 'fastq')
    inFq2 = SeqIO.parse(fq2, 'fastq')
    outputHandle1 = open(outFqPair1, 'w')
    outputHandle2 = open(outFqPair2, 'w')
    workingBarcode = ""
    workingSet1 = []
    workingSet2 = []
    for fqRec in inFq1:
        fqRec2 = inFq2.next()
        bc4fq1 = fqRec.description.split("###")[-2].strip()
        # print("Working barcode: {}. Current barcode: {}.".format(
        #       workingBarcode,barcodeRecord))
        # print("name of read with this barcode: {}".format(
        #       record.qname))
        # print("Working set: {}".format(workingSet1))
        a = "AAAAAAAAAA"
        t = "TTTTTTTTTT"
        g = "GGGGGGGGGG"
        c = "CCCCCCCCCC"
        if(a in bc4fq1 or t in bc4fq1 or c in bc4fq1 or g in bc4fq1):
            continue
        # if(int(fqRec.description.split('###')[-1].strip()) < 2):
        #    continue
        # Originally removing reads with family size <2, since one pair could
        # have more than the other, it's important that I keep these reads in
        # and filter them from the BAM file
        if(workingBarcode == ""):
            workingBarcode = bc4fq1
            workingSet1 = []
            workingSet1.append(fqRec)
            workingSet2 = []
            workingSet2.append(fqRec2)
            continue
        elif(workingBarcode == bc4fq1):
            workingSet1.append(fqRec)
            workingSet2.append(fqRec2)
            continue
        elif(workingBarcode != bc4fq1):
            mergedRecord1, success1 = compareFastqRecords(
                workingSet1, stringency=float(stringency))
            mergedRecord2, success2 = compareFastqRecords(
                workingSet2, stringency=float(stringency))
            if(success1):
                SeqIO.write(mergedRecord1, outputHandle1, "fastq")
            if(success2):
                SeqIO.write(mergedRecord2, outputHandle2, "fastq")
            workingSet1 = []
            workingSet1.append(fqRec)
            workingSet2 = []
            workingSet2.append(fqRec2)
            workingBarcode = bc4fq1
            continue
        else:
            raise RuntimeError(
                "No idea what's going on. This code should be unreachable")
    inFq1.close()
    inFq2.close()
    outputHandle1.close()
    outputHandle2.close()
    return outFqPair1, outFqPair2


def reverseComplement(fq, dest="default"):
    if(dest == "default"):
        temp = '.'.join(fq.split('.')[0:-1]) + "_rc.fastq"
        dest = temp.split('/')[-1]
    InFastq = SeqIO.parse(fq, "fastq")
    OutFastq = open(dest, "w", 0)
    for record in InFastq:
        rc_record = record.reverse_complement(id=record.id + "_rc")
        SeqIO.write(rc_record, OutFastq, "fastq")
    OutFastq.close()
    return dest


def singleFastqConsolidate(fq, outFq="default", stringency=0.9):
    if(outFq == "default"):
        outFq = '.'.join(fq.split('.')[0:-1]) + 'cons.fastq'
    inFq = SeqIO.parse(fq, 'fastq')
    outputHandle = open(outFq, 'w')
    workingBarcode = ""
    workingSet = []
    for fqRec in inFq:
        bc4fq = fqRec.description.split("###")[-2].strip()
        # print("Working barcode: {}. Current barcode: {}.".format(
        #   workingBarcode,barcodeRecord))
        # print("name of read with this barcode: {}".format(record.qname))
        # print("Working set: {}".format(workingSet))
        a = "AAAAAAAAAA"
        t = "TTTTTTTTTT"
        g = "GGGGGGGGGG"
        c = "CCCCCCCCCC"
        if(a in bc4fq or t in bc4fq or c in bc4fq or g in bc4fq):
            continue
        if(int(fqRec.description.split('###')[-1].strip()) < 2):
            continue
        if(workingBarcode == ""):
            workingBarcode = bc4fq
            workingSet = []
            workingSet.append(fqRec)
            continue
        elif(workingBarcode == bc4fq):
            workingSet.append(fqRec)
            continue
        elif(workingBarcode != bc4fq):
            mergedRecord, success = compareFastqRecords(
                workingSet, stringency=float(stringency))
            if(success):
                SeqIO.write(mergedRecord, outputHandle, "fastq")
            workingSet = []
            workingSet.append(fqRec)
            workingBarcode = bc4fq
            continue
        else:
            raise RuntimeError(
                "No idea what's going on. This code should be unreachable")
    inFq.close()
    outputHandle.close()
    return outFq


def TrimHoming(
        fq,
        homing,
        trimfq="default",
        bar_len=12,
        tags_file="default",
        trim_err="default",
        start_trim=1):
    pl("TrimHoming: \"{}\" from {}".format(homing, fq))
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    if(trim_err == "default"):
        temp = '.'.join(fq.split('.')[0:-1]) + '.err.fastq'
        trim_err = temp.split('/')[-1]
    if(trimfq == "default"):
        temp = '.'.join(fq.split('.')[0:-1]) + '.trim.fastq'
        trimfq = temp.split('/')[-1]
    if(tags_file == "default"):
        temp = '.'.join(fq.split('.')[0:-1]) + '.tags.fastq'
        tags_file = temp.split('/')[-1]
    tagsOpen = open(tags_file, "w", 0)
    trimOpen = open(trimfq, "w", 0)
    InFastq = SeqIO.parse(fq, "fastq")
    HomingLen = len(homing)
    TotalTrim = HomingLen + bar_len + start_trim
    pl("Homing Length is {}".format(HomingLen))
    for rec in InFastq:
        pre_tag = SeqRecord(
            Seq(str(rec.seq)[0:bar_len], "fastq"),
            id=rec.id, description=rec.description)
        pre_tag.letter_annotations['phred_quality'] = rec.letter_annotations[
            'phred_quality'][0:bar_len]
        '''
        if homing not in pre_tag.seq:
            pl("Homing sequence not in tag. Writing to error file.")
            SeqIO.write(rec,errOpen,"fastq")
            continue
        '''
        SeqIO.write(pre_tag, tagsOpen, "fastq")
        descString = rec.description + " ###" + str(rec.seq[0:bar_len])
        post_tag = SeqRecord(Seq(str(rec.seq)[TotalTrim:],
                                 "fastq"),
                             id=rec.id,
                             description=descString)
        post_tag.letter_annotations['phred_quality'] = rec.letter_annotations[
            'phred_quality'][TotalTrim:]
        SeqIO.write(post_tag, trimOpen, "fastq")
    tagsOpen.close()
    trimOpen.close()
    return(tags_file, trimfq)
