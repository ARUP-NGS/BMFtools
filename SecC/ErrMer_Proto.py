from __future__ import division
import numpy as np
from utilBMF.HTSUtils import pFastqProxy, pFastqFile

consInFq1 = pFastqFile("/home/brett/Projects/BMFTools_Devel/lamda_data/lamda-50di"
                    "v_S4_L001_R1_001.fastq.rescued.shaded.BS.fastq")
consInFq2 = pFastqFile("/home/brett/Projects/BMFTools_Devel/lamda_data/kmer_test/"
                    "")
