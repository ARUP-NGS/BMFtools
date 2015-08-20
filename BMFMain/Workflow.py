import subprocess
import time
import logging
from subprocess import check_call, CalledProcessError
try:
    from re2 import findall
except ImportError:
    print("re2 import failed - falling back to re.")
    from re import findall

import cython
import numpy as np

from MawCluster import BCBam, BCFastq
from MawCluster.BCFastq import TrimHomingPaired
from utilBMF import HTSUtils
from utilBMF.ErrorHandling import ThisIsMadness as Tim
from utilBMF.HTSUtils import printlog as pl, TrimExt, PipeAlignTag
from utilBMF.MPA import MPA2Bam
from utilBMF.QC import GetAllQCMetrics, GetFamSizeStats


@cython.locals(calcCoverage=cython.bint, coverageForAllRegions=cython.bint,
               addRG=cython.bint, rLen=int)
def pairedBamProc(consfq1, consfq2, opts="",
                  bamPrefix="default", ref="default", aligner="default",
                  bed="/yggdrasil/workspace/Barcode_Noah/cfDNA_targets.bed",
                  abrapath="default",
                  coverageForAllRegions=False,
                  calcCoverage=True,
                  bwapath="default",
                  picardpath="default",
                  addRG=False,
                  realigner="abra", gatkpath="default", dbsnp="default",
                  rLen=-1, intelDeflator="default",
                  kmers_precomputed=False,
                  sortMem="6G"):
    """
    Performs alignment and bam tagging of consolidated fastq files.
    """
    if(ref == "default"):
        raise ValueError("Reference index required!")
    if(bamPrefix == "default"):
        bamPrefix = '.'.join(consfq1.split('.')[0:-1])
    if(aligner == "default"):
        pl("No aligner set, defaulting to bwa mem.")
        aligner = "mem"
    if("gatk" in realigner.lower()):
        pl("Since realigning with GATK IndelRealigner, "
           "addRG is being set to True.")
        addRG = True
    if(aligner == "mem"):
        pl("Now aligning with: %s" % aligner)
        '''
        taggedBAM = BCBam.AlignAndTagMem(
            consfq1, consfq2, ref=ref, opts=opts)
        '''
        coorSortFull = PipeAlignTag(consfq1, consfq2, ref=ref,
                                    path=bwapath,
                                    coorsort=True, u=False,
                                    sortMem=sortMem, opts=opts)
    elif(aligner == "aln"):
        outBAMProperPair = HTSUtils.align_bwa_aln(consfq1, consfq2,
                                                  ref=ref, opts=opts,
                                                  addRG=addRG)
        pl("Now tagging BAM with custom SAM tags.")
        taggedBAM = BCBam.pairedBarcodeTagging(
            consfq1, consfq2, outBAMProperPair, realigner=realigner)
        coorSortFull = HTSUtils.CoorSortAndIndexBam(taggedBAM)
    else:
        raise ValueError("Sorry, only bwa is supported currently.")
    if(rLen < 0):
        pl("rLen < 0 in pairedBamproc. This typically means that this "
           "is not set.")
    if(picardpath == "default"):
        pl("Warning: path to picard jar not set. This isn't required for much"
           ", but in case something dies later, this could be responsible",
           level=logging.DEBUG)

    pl("Now realigning with: %s" % realigner)
    if("abra" in realigner.lower()):
        if(abrapath == "default"):
            raise Tim("abrapath must be set for abra to be the realigner")
        realignedFull = BCBam.AbraCadabra(coorSortFull, ref=ref, bed=bed,
                                          jar=abrapath, rLen=rLen,
                                          intelPath=intelDeflator,
                                          kmers_precomputed=kmers_precomputed)
        if("gatk" in realigner.lower()):
            pl("Realigning around known indels, too. "
               "Lots of steps, hard to know which ones matter.")
            realignedFull = BCBam.GATKIndelRealignment(coorSortFull, ref=ref,
                                                       bed=bed, gatk=gatkpath,
                                                       dbsnp=dbsnp)
    elif("gatk" in realigner.lower().split(",")):
        realignedFull = BCBam.GATKIndelRealignment(coorSortFull, ref=ref,
                                                   bed=bed, gatk=gatkpath,
                                                   dbsnp=dbsnp)
    else:
        pl("Skipping realignment.")
        realignedFull = coorSortFull
        # realignedFull = BCBam.AbraCadabra(taggedBAM, ref=ref, bed=bed)
    namesortedRealignedFull = HTSUtils.NameSort(realignedFull, uuid=True)

    # SVOutputFile = BCBam.CallTranslocations(SVBam, bedfile=bed)
    coorSorted = HTSUtils.CoorSortAndIndexBam(namesortedRealignedFull)
    check_call(["rm", namesortedRealignedFull])
    return coorSorted


def pairedFastqParallel(inFastq1, inFastq2, indexFq="default",
                        p3Seq="default", p5Seq="default",
                        overlapLen=6, sortMem="6G",
                        homing=None, bcLen=-1, head=0, rescue=False,
                        minFamRsq=10, mmRsq=1,
                        nbases=-1):
    pl("Beginning pairedFastqParallel for {}, {}".format(inFastq1, inFastq2))
    if(rescue):
        raise NotImplementedError(
            "pairedFastqParallel does not yet support"
            " rescue. Ask me for it!")
    else:
        SplitFqSets = BCFastq.PairedShadeSplitter(
            inFastq1, inFastq2, indexFq=indexFq,
            head=head, nbases=nbases)
    '''
    BSortFq1, BSortFq2 = BCFastq.BarcodeSortBoth(bcFastq1, bcFastq2,
                                                 sortMem=sortMem)
    BConsFastq1, BConsFastq2 = BCFastq.pairedFastqConsolidate(
        BSortFq1, BSortFq2)
    pl("Parameters for cutadapt are p3Seq={}, p5Seq={}".format(p3Seq, p5Seq))
    if(p3Seq != "default"):
        pl("Running cutadapt ...")
        BConsFastq1, BConsFastq2 = BCFastq.CutadaptPaired(
            BConsFastq1, BConsFastq2, overlapLen=overlapLen,
            p3Seq=p3Seq, p5Seq=p5Seq)
    else:
        pl("Skipping cutadapt ...")
    # check_call(["rm", BSortFq1, BSortFq2])
    famStats = GetFamSizeStats(
        BConsFastq1,
        outfile=TrimExt(inFastq1) + ".famstats.txt")
    '''
    return SplitFqSets


def pairedFastqShades(inFastq1, inFastq2, indexFq="default",
                      p3Seq="default", p5Seq="default",
                      overlapLen=6, sortMem="6G", inline_barcodes=False,
                      homing=None, bcLen=-1, head=0, rescue=False,
                      minFamRsq=10, mmRsq=1):
    pl("Beginning pairedFastqShades for {}, {}".format(inFastq1, inFastq2))
    if(inline_barcodes is False):
        if(rescue):
            bcFastq1, bcFastq2 = BCFastq.RescueShadingWrapper(
                inFastq1, inFastq2, indexFq=indexFq,
                minFam=minFamRsq, head=head, mm=mmRsq)
        else:
            bcFastq1, bcFastq2 = BCFastq.FastqPairedShading(inFastq1,
                                                            inFastq2,
                                                            indexFq=indexFq,
                                                            head=head)
    else:
        bcFastq1, bcFastq2 = TrimHomingPaired(inFastq1, inFastq2,
                                              homing=homing, bcLen=bcLen)
    BSortFq1, BSortFq2 = BCFastq.BarcodeSortBoth(bcFastq1, bcFastq2,
                                                 sortMem=sortMem)
    # check_call(["rm", bcFastq1, bcFastq2])
    BConsFastq1, BConsFastq2 = BCFastq.pairedFastqConsolidate(
        BSortFq1, BSortFq2)
    pl("Parameters for cutadapt are p3Seq={}, p5Seq={}".format(p3Seq, p5Seq))
    if(p3Seq != "default"):
        pl("Running cutadapt ...")
        BConsFastq1, BConsFastq2 = BCFastq.CutadaptPaired(
            BConsFastq1, BConsFastq2, overlapLen=overlapLen,
            p3Seq=p3Seq, p5Seq=p5Seq)
    else:
        pl("Skipping cutadapt ...")
    # check_call(["rm", BSortFq1, BSortFq2])
    famStats = GetFamSizeStats(
        BConsFastq1,
        outfile=TrimExt(inFastq1) + ".famstats.txt")
    return BConsFastq1, BConsFastq2


def singleFastqShades(inFastq, indexFq="default",
                      p3Seq="default", p5Seq="default",
                      overlapLen=6, sortMem="6G",
                      inline_barcodes=False, homing=None,
                      bcLen=-1, head=0):
    pl("Beginning singleFastqShades for {}".format(inFastq))
    if(inline_barcodes is False):
        bcFastq = BCFastq.FastqSingleShading(inFastq, indexFq=indexFq,
                                             head=head)
    else:
        bcFastq = BCFastq.TrimHomingSingle(inFastq, homing=homing, bcLen=bcLen)
    BSortFq = BCFastq.BarcodeSort(bcFastq, sortMem=sortMem)
    BConsFastq = BCFastq.singleFastqConsolidate(BSortFq)
    pl("Parameters for cutadapt are p3Seq={}, p5Seq".format(p3Seq, p5Seq))
    if(p3Seq != "default"):
        pl("Running cutadapt ...")
        BConsFastq = BCFastq.CutAdaptSingle(BConsFastq, overlapLen=overlapLen,
                                            p3Seq=p3Seq, p5Seq=p5Seq)
    else:
        pl("Skipping cutadapt")
    check_call(["rm", BSortFq])
    famStats = GetFamSizeStats(
        BConsFastq,
        outfile=TrimExt(inFastq) + ".famstats.txt")
    return BConsFastq
