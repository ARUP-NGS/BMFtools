import argparse

from utilBMF.HTSUtils import printlog as pl, IllegalArgumentError
from utilBMF import HTSUtils
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from MawCluster import BCBam
'''
Translocation Informatics by Fractionation, IE
XLIF: X-Location Informatics by Fractional Alignment
'''


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--paired',
        '-p',
        help=(
            "State whether the experiment is paired-end or single-end.\n"
            "True - paired. "
            "False - single. "),
        default="True")
    parser.add_argument(
        '-i',
        '--input-fq',
        nargs="+",
        required=True)
    parser.add_argument(
        '-n',
        '--nmer',
        help=("Select size for n-mer - must be a positive integer."
              "\nDefault: 30"),
        default="30")
    parser.add_argument(
        '--read-length',
        '-l',
        help=("Provide read length. Assumption of"
              "Illumina data and uniform length."
              "\nIf not set, will guess from the data."),
        default="default"
    )
    parser.add_argument(
        '--name_tag',
        '-t',
        help="Optional tag for output read ids.",
        default="default")
    parser.add_argument("--barcodeDescTag",
                        "-b",
                        help=("Barcode Description Tag."
                              "\nIf this argument is not default,"
                              " then all information contained in the"
                              " read description tagged with these characters"
                              " will be kept in the split reads."
                              "\nAdditionally, a \"PartN\" string will be"
                              " appended to the description to avoid confusion"
                              " with the original fastqs")
                        )
    parser.add_argument("--conf",
                        "-c",
                        help=("Configuration file path containing keys and"
                              "values separated by \"=\"."
                              ),
                        )
    parser.add_argument("--ref",
                        '-r',
                        help="Full path to reference file/index base.",
                        default="default")
    parser.add_argument("--alignment_opts",
                        "-a",
                        help="Additional alignment options for alignment.",
                        default="")
    parser.add_argument("--outPrefix",
                        "-o",
                        help="Prefix for out files. Defaults based on input",
                        default="default")
    parser.add_argument("--picardDir",
                        help="Directory for Picard Jar files.",
                        default="default")
    args = parser.parse_args()
    fq = args.input_fq
    paired = True
    '''
    if(args.ref == "default" and args.conf == "default"):
        raise IllegalArgumentError("You must provide either a refe")
    '''
    ref = args.ref
    if(args.outPrefix != "default"):
        outPrefix = args.outPrefix
    else:
        outPrefix = '.'.join(fq[0].split('.'))
    if(args.paired.lower() != "true"):
        paired = False
        pl("Single-end selected.")
        if(len(fq) != 1):
            raise IllegalArgumentError(
                "Single-end requires exactly one fastq file.")
    else:
        pl("Paired-end selected.")
        if(len(fq) != 2):
            raise IllegalArgumentError(
                "Paired-end requires exactly two fastq files.")
    if(paired):
        pl("Opening file buffers")
        if(args.name_tag == "default"):
            outFqs = pairedSplitToNmers(fq[0],
                                        fq[1],
                                        read_length=args.read_length,
                                        n=args.nmer)
        else:
            outFqs = pairedSplitToNmers(fq[0],
                                        fq[1],
                                        read_length=args.read_length,
                                        n=args.nmer,
                                        name_tag=args.name_tag)
    else:
        pl("Opening file buffers.")
        if(args.name_tag == "default"):
            outFqs = SplitFastqToNmerReads(fq[0],
                                           read_length=args.read_length,
                                           n=args.nmer)
        else:
            outFqs = SplitFastqToNmerReads(fq[0],
                                           read_length=args.read_length,
                                           n=args.nmer,
                                           name_tag=args.name_tag)
        pl("File buffers closed.")
    outBAMs = []
    if(len(outFqs[1]) == 1):
        for i in outFqs:
            HTSUtils.align_bwa_se(i, ref, args.alignment_opts,
                                  outsam='.'.join(i.split('.')[0:1]) + '.sam')
            BCBam.Sam2Bam('.'.join(
                i.split('.')[0:-1]) + '.sam', '.'.join(
                    i.split('.')[0:-1]) + '.bam')
            outBAMs.append('.'.join(i.split('.')[0:1]) + '.bam')
    else:
        for i in outFqs:
            HTSUtils.align_bwa(i[0], i[-1], ref, args.alignment_opts,
                               outsam='.'.join(i[0].split('.')[0:-1]) + '.sam')
            BCBam.Sam2Bam('.'.join(
                i[0].split('.')[0:-1]) + '.sam', '.'.join(
                    i[0].split('.')[0:-1]) + '.bam')
            outBAMs.append('.'.join(i[0].split('.')[0:-1]) + '.bam')
    if(args.PicardDir != "default"):
        outMerged = HTSUtils.mergeBam(outBAMs,
                                      MergeJar=(args.PicardDir +
                                                "MergeSamFiles.jar"))
    return outMerged


def SplitFastqToNmerReads(inFastq, read_length="default", n="30",
                          name_tag="Split_Pt#",
                          outPrefix="default",
                          barcodeDescTag="default"):
    """Takes as input a fastq file, a number n for the size
    of the reads to be cut down to, and an optional
    name_tag, which is added to the end of the read id, followed by
    the piece number provided.
    """
    pl("Beginning splitting process for fastq. If data is paired-end,"
       " then this will be called twice.")
    if(outPrefix == "default"):
        outPrefix = '.'.join(inFastq.split(".")[0:-1])
    if(read_length == "default"):
        pl("Read length not specified, attempting to guess.")
        length = len(SeqIO.parse(inFastq, "fastq").next().seq)
        pl("Read length inferred as being " + str(length))
    else:
        read_length = int(read_length)
    try:
        n = int(n)
    except(ValueError):
        pl(("Given string for nmer size could not be converted to int."
           " Resetting to default value of 30"))
        n = 30
    numOutputs = read_length/n
    outBuffers = []
    outFiles = [
        outPrefix + ".N" + str(n) + ".Pt" + str(i) + ".fastq" for i in range(
            numOutputs)]
    for i in outFiles:
        outBuffers.append(open(i, "w"))
    for seqRec in SeqIO.parse(inFastq, "fastq"):
        count = 1
        for buffer in outBuffers:
            newDesc = ""
            if(barcodeDescTag == "default"):
                newDesc = seqRec.id + name_tag + str(
                    count) + seqRec.description[
                    len(seqRec.name):]
            else:
                oldDescList = [i.strip() for i in seqRec.description.split(
                    barcodeDescTag)]
                if(len(oldDescList) != 4):
                    pl("Didn't find provided description split tag. Ignoring.")
                    newDesc = seqRec.id + name_tag + str(
                        count) + seqRec.description[
                        len(seqRec.name):]
                else:
                    # Then add a "PartN", where N is the number, to the barcode
                    # to eliminate risks of merging families from different
                    # portions of the read
                    newDesc = ((" " + barcodeDescTag).join(oldDescList[0:-2]) +
                               " " + barcodeDescTag + oldDescList[-2] + "Part"
                               + str(count) + " " + barcodeDescTag +
                               oldDescList[-1])
            if(count != len(outBuffers)):
                tempSeq = SeqRecord(seq=seqRec.seq[(count-1)*n:count*n],
                                    id=seqRec.id + name_tag + str(count),
                                    name=seqRec.name + name_tag + str(count),
                                    description=seqRec.id + name_tag + str(
                                    count) + seqRec.description[
                                    len(seqRec.name):])
                tempSeq.letter_annotations[
                    'phred_quality'] = seqRec.letter_annotations[
                        'phred_quality'][count*n-n:count*n]
            else:
                tempSeq = SeqRecord(seq=seqRec.seq[count*n-n:],
                                    id=seqRec.id + name_tag + str(count),
                                    name=seqRec.name + name_tag + str(count),
                                    description=seqRec.id + name_tag + str(
                                    count) + seqRec.description[
                                    len(seqRec.name):])
                tempSeq.letter_annotations[
                    'phred_quality'] = seqRec.letter_annotations[
                        'phred_quality'][count*n-n:]
            SeqIO.write(tempSeq, buffer, "fastq")
            count += 1
    for buffer in outBuffers:
        buffer.close()
    return outFiles


def pairedSplitToNmers(inFastq1, inFastq2, read_length="default", n="30",
                       name_tag="Split_Pt#",
                       outPrefix="default",
                       barcodeDescTag="default"):
    outFilesR1 = SplitFastqToNmerReads(inFastq1,
                                       read_length,
                                       n,
                                       name_tag,
                                       outPrefix,
                                       barcodeDescTag)
    outFilesR2 = SplitFastqToNmerReads(inFastq2,
                                       read_length,
                                       n,
                                       name_tag,
                                       outPrefix,
                                       barcodeDescTag)
    return zip(outFilesR1, outFilesR2)

if(__name__ == "__main__"):
    main()
