from Bio import SeqIO
from HTSUtils import printlog as pl

"""
Contains various utilities for working with barcoded fastq files.

TODO: change steps relying on pass/fail for homing sequence to check only
for presence or absence of "pass" in the .lower() of that field,
as a pass/fail will be used for the i7 indexing protocol to take
its place, generalizing the code.

TODO: Try using a preallocated numpy multidimensional array to speed up
inexact matching. compareFastqRecordsInexactNumpy

"""


def BarcodeSort(inFastq, outFastq="default"):
    pl("Sorting {} by barcode sequence.".format(inFastq))
    from subprocess import call
    if(outFastq == "default"):
        outFastq = '.'.join(inFastq.split('.')[0:-1]) + '.BS.fastq'
    BS1 = ("cat " + inFastq + " | paste - - - - | grep 'HomingPass' | sed "
           "'s:###:###\t:g' |  awk 'BEGIN {{FS=OFS=\"\t\"}};{{print $3,$0}}' |"
           " sort -k1,1 | awk 'BEGIN {{OFS=FS=\"\t\"}};"
           "{{print $2,$3,$4,$5,$6,$7,$8,$9,$10}}' | "
           "sed 's:###\t:###:g' | sed 's:\t$::g' "
           "| sed 's:\t$::g' | tr '\t' '\n' > ")
    BSstring = BS1 + outFastq
    call(BSstring, shell=True)
    pl("Command: {}".format(BSstring.replace(
        "\t", "\\t").replace("\n", "\\n")))
    return outFastq


def compareFastqRecords(R, stringency=0.9, hybrid=False):
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
    elif(hybrid is True):
        return compareFastqRecordsInexact(R)
    consolidatedRecord = SeqRecord(seq=finalSeq, id=R[0].id,
                                   letter_annotations=R[0].letter_annotations,
                                   name=R[0].name,
                                   description=R[0].description)
    return consolidatedRecord, Success


def compareFastqRecordsInexact(R):
    import math
    import numpy as np
    """Uses the phred scores of base qualities to """
    # pl("Doing inexact fastq merging for family.")
    from Bio.SeqRecord import SeqRecord
    seqs = [str(record.seq) for record in R]
    quals = [record.letter_annotations['phred_quality'] for record in R]
    Success = True
    newQuals = []
    newSeq = ""
    count = 0
    for i in range(len(seqs[0])):
        candidates = zip([s[i] for s in seqs], [q[i] for q in quals])
        candidates.sort(key=lambda tup: tup[0])
        # print(repr(candidates) + " is candidates repr.")
        canBases = list(set([c[0] for c in candidates]))
        canBases.sort()
        # print(repr(canBases) + " is canBases repr.")
        pIncorr = []  # 4 x 100 array, tallying p-scores ?
        for base in canBases:
            bSupport = [c for c in candidates if c[0] == base]
            pIncorr.append(prod([math.pow(10, -c[1] / 10.) for c in bSupport]))
        if(len(pIncorr) > 1):
            winners = zip(canBases, pIncorr)
            winners.sort(key=lambda tup: tup[1])
            winner = winners[0]
            del winners
        else:
            winner = zip(canBases, pIncorr)[0]
        newSeq += (winner[0])
        newQuals.append(
            int(-10 * np.log10(winner[1] * winner[1] / np.prod(pIncorr))))
        count += 1
        # if(count < 20):
        #     print(str(newQuals[-1]) + " = phred score for this base.")
    newQuals = [i if i <= 93 else 93 for i in newQuals]
    consolidatedRecord = SeqRecord(
        seq=newSeq,
        id=R[0].id,
        name=R[0].name,
        description=R[0].name)
    consolidatedRecord.letter_annotations['phred_quality'] = newQuals
    return consolidatedRecord, Success


def compareFastqRecordsInexactNumpy(R):
    import numpy as np
    from Bio.SeqRecord import SeqRecord
    letterNumDict = {}
    letterNumDict['A'] = 0
    letterNumDict['C'] = 1
    letterNumDict['G'] = 2
    letterNumDict['T'] = 3
    letterNumDict[0] = 'A'
    letterNumDict[1] = 'C'
    letterNumDict[2] = 'G'
    letterNumDict[3] = 'T'

    def dAccess(x):
        return letterNumDict[x]
    dAccess = np.vectorize(dAccess)
    seqs = [str(record.seq) for record in R]
    stackArrays = tuple([np.char.array(s, itemsize=1) for s in seqs])
    seqArray = np.vstack(stackArrays)
    # print(repr(seqArray))
    quals = [record.letter_annotations['phred_quality'] for record in R]
    qualProb = np.power(10, np.multiply(-1/10., quals))
    qualA = qualProb.copy()
    qualC = qualProb.copy()
    qualG = qualProb.copy()
    qualT = qualProb.copy()
    qualA[seqArray != "A"] = 1.0
    qualAProd = np.prod(qualA, 0)
    qualC[seqArray != "C"] = 1.0
    qualCProd = np.prod(qualC, 0)
    qualG[seqArray != "G"] = 1.0
    qualGProd = np.prod(qualG, 0)
    qualT[seqArray != "T"] = 1.0
    qualTProd = np.prod(qualT, 0)
    qualAllProd = np.vstack([qualAProd, qualCProd, qualGProd, qualTProd])
    '''
    v2 - maybe faster?
    qualA[seqArray != "A"] = 1.0
    qualC[seqArray != "C"] = 1.0
    qualG[seqArray != "G"] = 1.0
    qualT[seqArray != "T"] = 1.0
    qualAllProd = np.vstack([
        np.prod(qualA, 0), np.prod(
        qualC, 0), np.prod(qualG, 0), np.prod(qualT, 0)])

    '''
    Success = True
    newSeq = "".join(
        np.apply_along_axis(dAccess, 0, np.argmin(qualAllProd, 0)))
    # Calculates new probability that the best
    # candidate is right and all the others are wrong
    # based on base qualities.
    # NOT WORKING
    # phredQuals = np.multiply(
    #    -10, np.log10(np.multiply(np.power(np.prod(
    #         qualAllProd, 0), -1), np.power(
    #         np.amin(qualAllProd, 0), 2)))).astype(int)
    # omgz = np.multiply(np.power(np.prod(
    # qualAllProd, 0), -1), np.power(np.amin(qualAllProd, 0), 2))
    # if((phredQuals < 0).any()):
    #    print(repr(qualAllProd))
    #    print(repr(np.prod(qualAllProd, 0)))
    #    print(repr(np.amin(qualAllProd, 0)))
    #    raise ValueError("Shouldn't be negative.")
    phredQuals = np.multiply(
        -10, np.log10(np.amin(qualAllProd, 0))).astype(int)
    consolidatedRecord = SeqRecord(
        seq=newSeq,
        id=R[0].id,
        name=R[0].name,
        description=R[0].description)
    consolidatedRecord.letter_annotations[
        'phred_quality'] = [i if i <= 93 else 93 for i in phredQuals.tolist()]
    return consolidatedRecord, Success


def FastqPairedShading(fq1,
                       fq2,
                       indexfq="default",
                       outfq1="default",
                       outfq2="default",
                       gzip=True):
    pl("Now beginning fastq marking: Pass/Fail and Barcode")
    if(indexfq == "default"):
        raise ValueError("For an i5/i7 index ")
    if(outfq1 == "default"):
        outfq1 = '.'.join(fq1.split('.')[0:-1]) + '.shaded.fastq'
    if(outfq2 == "default"):
        outfq2 = '.'.join(fq1.split('.')[0:-1]) + '.shaded.fastq'
    inFq1 = SeqIO.parse(fq1, "fastq")
    inFq2 = SeqIO.parse(fq2, "fastq")
    outFqHandle1 = open(outfq1, "w")
    outFqHandle2 = open(outfq2, "w")
    inIndex = SeqIO.parse(indexfq)
    for read1 in inFq1:
        read2 = inFq2.next()
        indexRead = inIndex.next()
        if("N" in indexRead.seq):
            read1.description += " ###IndexFail ###" + indexRead.seq
            read2.description += " ###IndexFail ###" + indexRead.seq
            SeqIO.write(read1, outFqHandle1, "fastq")
            SeqIO.write(read2, outFqHandle2, "fastq")
        else:
            read1.description += " ###IndexPass ###" + indexRead.seq
            read2.description += " ###IndexPass ###" + indexRead.seq
            SeqIO.write(read1, outFqHandle1, "fastq")
            SeqIO.write(read2, outFqHandle2, "fastq")
    outFqHandle1.close()
    outFqHandle2.close()
    if(gzip is True):
        from subprocess import call
        call(['gzip', fq1], shell=False)
        call(['gzip', fq2], shell=False)
    return outfq1, outfq2


def FastqSingleShading(fq,
                       indexfq="default",
                       outfq="default",
                       gzip=True):
    pl("Now beginning fastq marking: Pass/Fail and Barcode")
    if(indexfq == "default"):
        raise ValueError("For an i5/i7 index ")
    if(outfq == "default"):
        outfq = '.'.join(fq.split('.')[0:-1]) + '.shaded.fastq'
    inFq1 = SeqIO.parse(fq, "fastq")
    outFqHandle1 = open(outfq, "w")
    inIndex = SeqIO.parse(indexfq)
    for read1 in inFq1:
        indexRead = inIndex.next()
        if("N" in indexRead.seq):
            read1.description += " ###IndexFail ###" + indexRead.seq
            SeqIO.write(read1, outFqHandle1, "fastq")
        else:
            read1.description += " ###IndexPass ###" + indexRead.seq
            SeqIO.write(read1, outFqHandle1, "fastq")
    outFqHandle1.close()
    if(gzip is True):
        from subprocess import call
        call(['gzip', fq], shell=False)
    return


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
    command_str = ('cat ' + in_fastq + ' | paste - - - - | '
                   'sort -k1,1 -t " " | tr "\t" "\n"')
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
    CommandStr = ("cat " + fq + " | paste - - - - | grep '" +
                  string + "' | tr '\t' '\n' > " + matchFile)
    call(CommandStr, shell=True)
    CommandStr2 = ("cat {} | paste - - - - | grep ".format(fq) +
                   "-v '{}' | tr '\t' '\n' > {}".format(
        string,
        missFile))
    call(CommandStr2, shell=True)
    return(CommandStr, CommandStr2, matchFile, missFile)


def fastx_trim(infq, outfq, n):
    pl("Now beginning fastx_trimmer.")
    import subprocess
    command_str = ['fastx_trimmer', '-l', str(n), '-i', infq, '-o', outfq]
    pl(command_str)
    subprocess.call(command_str)
    return(command_str)


def getProperPairs(infq1, infq2, shared="default", outfq1="default",
                   outfq2="default", outfqSingle="default"):
    """Assuming that the information here is sorted by barcode..."""
    pl("Now beginning getProperPairs.")
    from Bio import SeqIO
    if(shared == "default"):
        raise ValueError(
            "An index for the barcodes shared must be present.")
    if(outfq1 == "default"):
        outfq1 = '.'.join(infq1.split('.')[0:-1]) + '.proper.fastq'
    if(outfq2 == "default"):
        outfq2 = '.'.join(infq2.split('.')[0:-1]) + '.proper.fastq'
    if(outfqSingle == "default"):
        outfqSingle = '.'.join(infq1.split('.')[0:-1]) + '.solo.fastq'
    infq1Handle = SeqIO.parse(infq1, "fastq")
    infq2Handle = SeqIO.parse(infq2, "fastq")
    outfq1Handle = open(outfq1, 'w')
    outfq2Handle = open(outfq2, 'w')
    outfqSingleHandle = open(outfqSingle, 'w')

    f = open(shared, "r")
    BarDict = {}
    for line in f.readlines():
        BarDict[line.strip()] = ""
    f.close()
    for read in infq2Handle:
        try:
            BarDict[read.description.split('###')[-2].strip()]
            # Testing if it's in the dictionary, no need to assign
            # This is simply checking whether or not
            # it has a mate with the same barcode sequence.
            SeqIO.write(read, outfq2Handle, "fastq")
        except KeyError:
            SeqIO.write(read, outfqSingleHandle, "fastq")

    for read in infq1Handle:
        try:
            BarDict[read.description.split('###')[-2].strip()]
            SeqIO.write(read, outfq1Handle, "fastq")
        except KeyError:
            SeqIO.write(read, outfqSingleHandle, "fastq")
    return outfq1, outfq2, outfqSingle


def getSharedBC(barIdx1, barIdx2, shared="default"):
    from subprocess import call
    if(shared == "default"):
        shared = '.'.join(barIdx1.split('.')[0:-1]) + '.sharedBC'
    Str = "cat {} {} | awk '{{print $2}}' | sort |".format(barIdx1, barIdx2)
    Str += " uniq -c | grep -v '1' | awk '{{print $NF}}' > {}".format(shared)
    pl("Command for shared barcodes is: {}".format(Str))
    call(Str, shell=True)
    return shared


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


def GenerateShadesIndex(indexFastq, index_file="default"):
    from subprocess import call
    if(index_file == "default"):
        index_file = '.'.join(indexFastq.split('.')[0:-1]) + ".barIdx"
    commandStr = ("cat {} | paste - - - - | cut -f2 | ".format(indexFastq) +
                  "sort | uniq -c | awk 'BEGIN {{OFS=\"\t\"}};{{print $1"
                  ",$2}}' > {}".format(index_file))
    call(commandStr, shell=True)
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


def halveFqRecords(fq, outfq1="default", outfq2="default", minLength=40):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    if(outfq1 == "default"):
        outfq1 = '.'.join(fq.split('.')[0:-1]) + '_R1.fastq'
    if(outfq2 == "default"):
        outfq2 = '.'.join(fq.split('.')[0:-1]) + '_R2.fastq'
    infq = SeqIO.parse(fq, "fastq")
    out1 = open(outfq1, "w")
    out2 = open(outfq2, "w")
    for read in infq:
        if(len(read.seq) < minLength):
            continue
        halfLen = len(read.seq) / 2
        #  read2.id = read1.id
        #  read2.description = ' '.join(
        #      read1.description.split()[1:]).replace('1:N', '2:N')
        outread1 = SeqRecord(Seq(str(read.seq[0:halfLen]), "fastq"),
                             id=read.id,
                             description=read.description + ' 1:N:0:1')
        outread1.letter_annotations[
            'phred_quality'] = read.letter_annotations[
            'phred_quality'][0:halfLen]
        SeqIO.write(outread1, out1, "fastq")
        outread2 = SeqRecord(Seq(str(read.seq[halfLen:]), "fastq"),
                             id=read.id,
                             description=read.description + ' 2:N:0:1')
        outread2.letter_annotations[
            'phred_quality'] = read.letter_annotations[
            'phred_quality'][halfLen:]
        SeqIO.write(outread2, out2, "fastq")
    out1.close()
    out2.close()
    return


def GetFamilySizePaired(
        trimfq1,
        trimfq2,
        BarcodeIndex,
        outfq1="default",
        outfq2="default",
        singlefq1="default",
        singlefq2="default"):
    pl("Running GetFamilySizeSingle for {}, {}.".format(trimfq1, trimfq2))
    infq1 = SeqIO.parse(trimfq1, "fastq")
    infq2 = SeqIO.parse(trimfq2, "fastq")
    if(outfq1 == "default"):
        outfq1 = '.'.join(trimfq1.split('.')[0:-1]) + ".fam.fastq"
    if(outfq2 == "default"):
        outfq2 = '.'.join(trimfq2.split('.')[0:-1]) + ".fam.fastq"
    outfqBuffer1 = open(outfq1, "w")
    outfqBuffer2 = open(outfq2, "w")
    if(singlefq1 == "default"):
        singlefq1 = '.'.join(
            trimfq1.split('.')[
                0:-1]) + ".lonely.hearts.club.band.fastq"
    if(singlefq2 == "default"):
        singlefq2 = '.'.join(
            trimfq2.split('.')[
                0:-1]) + ".lonely.hearts.club.band.fastq"
    singlefqBuffer1 = open(singlefq1, "w")
    singlefqBuffer2 = open(singlefq2, "w")
    TotalReads, ReadsWithFamilies = 0, 0
    index = open(BarcodeIndex, "r")
    dictEntries = [line.split() for line in index]
    BarDict = {}
    for entry in dictEntries:
        BarDict[entry[1]] = entry[0]
    for read in infq1:
        read2 = infq2.next()
        # index.seek(0)
        TotalReads += 1
        # print("description is {}".format(read.description))
        # print("-1 index of description is {}".format(
        # read.description.split("###")[-1].strip()))
        # print("-2 index of description is {}".format(
        # read.description.split("###")[-2].strip()))
        readTag1 = read.description.split("###")[-1].strip()
        newRead1 = read
        newRead2 = read2
        # print("readTag is _" + readTag + "_")
        try:
            famSize = BarDict[readTag1]
        except KeyError:
            famSize = 0
        newRead1.description = read.description + " ###" + str(famSize)
        newRead2.description = read.description + " ###" + str(famSize)
        # print("famSize = _{}_".format(str(famSize)))
        # print("The value of this comparison to 1 is {}".format(
        # str(famSize=="1")))
        if(str(famSize) == "0"):
            continue
        if(str(famSize) == "1" or str(famSize) == "2"):
            #  print("Hey, I found a singleton")
            SeqIO.write(newRead1, singlefqBuffer1, "fastq")
            SeqIO.write(newRead2, singlefqBuffer2, "fastq")
        else:
            #  print("Hey, I found a read with a family!")
            ReadsWithFamilies += 1
            SeqIO.write(newRead1, outfqBuffer1, "fastq")
            SeqIO.write(newRead2, outfqBuffer2, "fastq")
    outfqBuffer1.close()
    outfqBuffer2.close()
    singlefqBuffer1.close()
    singlefqBuffer2.close()
    return ([outfq1, outfq2],
            [singlefq1, singlefq2],
            TotalReads,
            ReadsWithFamilies)


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
            #  print("Hey, I found a singleton")
            SeqIO.write(newRead, singlefqBuffer, "fastq")
        else:
            #  print("Hey, I found a read with a family!")
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


def PairFastqBarcodeIndex(taggedFile1, taggedFile2, index_file="default"):
    pl("Now beginning GenerateFullFastqBarcodeIndex for {} and {}.".format(
        taggedFile1, taggedFile2))
    from subprocess import call
    if(index_file == "default"):
        index_file = '.'.join(taggedFile1.split('.')[0:-1]) + ".barIdx"
    cmd = "cat {} {} | sed 's: ###:\t:g' | paste - - - - | awk ".format(
        taggedFile1, taggedFile2)
    cmd += "'BEGIN {{FS=\"\t\"}};{{print $3}}' | sort | uniq -c | awk 'BEGIN "
    cmd += "{{OFS=\"\t\"}};{{print $1,$2}}' | sort -k1,1n > {}".format(
        index_file)
    pl("CommandStr = {}".format(cmd.replace("\t", "\\t")))
    call(cmd, shell=True)
    return index_file


def pairedFastqConsolidate(fq1, fq2, outFqPair1="default",
                           outFqPair2="default",
                           stringency=0.9, inexact=False,
                           numpy=False):
    if(numpy is True and inexact is True):
        raise ValueError("You must choose either numpy or inexact, not both.")
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
        # and filter them from the BAM filea
        numFams = 0
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
            # pl("About to collapse a family.")
            if(inexact is False and numpy is False):
                mergedRecord1, success1 = compareFastqRecords(
                    workingSet1, stringency=float(stringency))
                mergedRecord2, success2 = compareFastqRecords(
                    workingSet2, stringency=float(stringency))
            elif(inexact is True):
                mergedRecord1, success1 = compareFastqRecords(
                    workingSet1, hybrid=True)
                mergedRecord2, success2 = compareFastqRecords(
                    workingSet2, hybrid=True)
            elif(numpy is True):
                mergedRecord1, success1 = compareFastqRecordsInexactNumpy(
                    workingSet1)
                mergedRecord2, success2 = compareFastqRecordsInexactNumpy(
                    workingSet2)
            if(success1):
                SeqIO.write(mergedRecord1, outputHandle1, "fastq")
            if(success2):
                SeqIO.write(mergedRecord2, outputHandle2, "fastq")
            workingSet1 = []
            workingSet1.append(fqRec)
            workingSet2 = []
            workingSet2.append(fqRec2)
            workingBarcode = bc4fq1
            numFams += 1
            continue
        else:
            raise RuntimeError(
                "No idea what's going on. This code should be unreachable")
    inFq1.close()
    inFq2.close()
    outputHandle1.close()
    outputHandle2.close()
    pl("Total families merged together: " + str(numFams))
    return outFqPair1, outFqPair2


def prod(iterable):
    import operator
    return reduce(operator.mul, iterable, 1)


def renameReads(fq1, fq2, outfq1="default", outfq2="default"):
    #  Requires barcode-sorted, consolidated fastq files,
    #  filtered to only the shared barcode families
    if(outfq1 == "default"):
        outfq1 = fq1.split('.')[0] + '.cons.R1fastq'
    if(outfq2 == "default"):
        outfq2 = fq2.split('.')[0] + '.cons.R2fastq'
    from Bio import SeqIO
    infq1 = SeqIO.parse(fq1, "fastq")
    infq2 = SeqIO.parse(fq2, "fastq")
    outfqhandle1 = open(outfq1, "w")
    outfqhandle2 = open(outfq2, "w")
    for read1 in infq1:
        read2 = infq2.next()
        read2.id = read1.id
        read2.description = ' '.join(
            read1.description.split()[1:]).replace('1:N', '2:N')
        SeqIO.write(read1, outfqhandle1, "fastq")
        SeqIO.write(read2, outfqhandle2, "fastq")
    outfqhandle1.close()
    outfqhandle2.close()
    infq1.close()
    infq2.close()
    return outfq1, outfq2


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


def ShadesRescuePaired(singlefq1,
                       singlefq2,
                       appendFq1="default",
                       appendFq2="default",
                       index="default"):
    if(appendFq1 == "default" or appendFq2 == "default"):
        raise ValueError("Fastq for appending rescued reads required.")
    if(index == "default"):
        raise ValueError("Index for initial family sizes required.")
    ind = open(index, "r")
    dictEntries = [line.split() for line in ind]
    BarDict = {}
    for entry in dictEntries:
        BarDict[entry[1]] = entry[0]
    a = BarDict.keys()
    pass
    return


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
