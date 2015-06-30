import sys
from copy import copy as ccopy
from utilBMF.ErrorHandling import IllegalArgumentError

"""
This module contains default values for each subcommand of bmftools.
When adding additional subcommands or extending existing ones, there is a
simple protocol to follow.
First, write argparse normally, except that no defaults should be set.
This causes the Namespace to populate None objects, which this defaultConfig
overrides.
Second, for attributes for which you want a default behavior besides None set,
set that value here, replacing any '-' characters with '_', and make sure it
is clearly typed. (See below)
Third, when overriding the default config with a run config file,
make sure it follows "new config" style, which is as follows.

key|value|typechar#Comment field
e.g.:
abrapath|default|s
maxAF|0.1|f

typechar is 's' for string, 'b' for bool, 'i' for int, 'f' for float.

In addition, you must add a type= argument for each argument that needs to be
parsed in as something other than a string.

"""


defaultConfig = {
                 "abrapath": "default",
                 "aligner": "mem",
                 "analysisTag": "default",
                 "bcLen": -1,
                 "bed": "default",
                 "bed_buffer": 0,
                 "bwapath": "bwa",
                 "check_both": True,
                 "compression": "bb",
                 "coorsort": False,
                 "experiment": "",
                 "file_prefix": "default",
                 "gatkpath": "default",
                 "genome_size": 3.2e9,
                 "head": 2,
                 "homing": "CATG",
                 "intelDeflator": "default",
                 "insert_distance": 35,
                 "is_slave": False,
                 "keepConsensus": False,
                 "ligation_efficiency": 0.5,
                 "logfile": "default",
                 "mapped_fraction": 0.83,
                 "maxAF": 0.1,
                 "MaxPValue": 1e-15,
                 "mean_aligned_fraction": 0.9,
                 "minBQ": 0,
                 "minClustDepth": 10,
                 "minCov": 5,
                 "minFA": 3,
                 "minFM": 3,
                 "minFam": 10,
                 "minFracAgreed": 0.75,
                 "minMQ": 0,
                 "minPileupLen": 10,
                 "mismatches": 0,
                 "mm": 1,
                 "on_target": 0.25,
                 "opts": "",
                 "outBAM": None,
                 "outfile": "default",
                 "outfile_handle": sys.stdout,
                 "outTsv": "default",
                 "outVCF": "default",
                 "parallel": False,
                 "p3Seq": "default",
                 "p5Seq": "default",
                 "padding": 0,
                 "padding_distance": 120,
                 "qc_fail": 0.1,
                 "readLength": -1,
                 "realigner": "abra",
                 "ref": "default",
                 "rescue": False,
                 "review_dir": "default",
                 "single_end": False,
                 "sortMem": "6G",
                 "threads": 4,
                 "uncompressed_bam": False
                 }
