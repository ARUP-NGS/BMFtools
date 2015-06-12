import sys

"""
This module contains default values for each subcommand of bmftools.
When adding additional subcommands or extending existing ones, there is a simple
protocol to follow.
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

"""


defaultConfig = {
                 "abrapath": "default",
                 "aligner": "mem",
                 "analysisTag": "default",
                 "bed": "/dev/random",
                 "bed_buffer": 0,
                 "check_both": True,
                 "gatkpath": "default",
                 "genome_size": 3.2e9,
                 "head": 2,
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
                 "minFA": 3,
                 "minFam": 10,
                 "min_frac_agreed": 0.75,
                 "minFracAgreed": 0.75,
                 "minMQ": 0,
                 "minPileupLen": 10,
                 "mismatches": 0,
                 "mm": 1,
                 "on_target": 0.25,
                 "outfile": "default",
                 "outfile_handle": sys.stdout,
                 "outTsv": "default",
                 "outVCF": "default",
                 "p3Seq": "default",
                 "p5Seq": "default",
                 "padding": 0,
                 "padding_distance": 120,
                 "qc_fail": 0.1,
                 "readlength": -1,
                 "realigner": "abra",
                 "ref": "default",
                 "rescue": False,
                 "threads": 4
                 }