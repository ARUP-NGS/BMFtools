import sys
from copy import copy as ccopy
from utilBMF.ErrorHandling import IllegalArgumentError, ImproperArgumentError

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

For cases where it is a more complex type, such as a file handle, use 'o'.
These types cannot be instantiated properly from command line parsing.

In addition, you must add a type= argument for each argument that needs to be
parsed in as something other than a string.

"""

'''
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
                 "create_log": False,
                 "experiment": "",
                 "file_prefix": "default",
                 "gatkpath": "default",
                 "genome_size": 3.2e9,
                 "head": 4,
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
                 "num_nucs": None,
                 "on_target": 0.25,
                 "opts": "",
                 "outBAM": None,
                 "outfile": "default",
                 "outfile_handle": sys.stdout,
                 "outTsv": "default",
                 "outVCF": "default",
                 "overlapLen": 6,
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
                 "slave_sort_mem": "768M",
                 "single_end": False,
                 "sortMem": "6G",
                 "threads": 4,
                 "uncompressed_bam": False
                 }
'''

def GetArgumentType(cystr key):
    try:
        return DefaultConfig[key][1]
    except KeyError:
        return 's'
        '''
        raise ImproperArgumentError("Key %s not found in _bmftools" % key +
                                    "_helper. Invalid argument!")
        '''


DefaultConfig = {
                 "abrapath": ("default", 's'),
                 "aligner": ("mem", 's'),
                 "analysisTag": ("default", 's'),
                 "bcLen": (-1, 'i'),
                 "bed": ("default", 's'),
                 "bed_buffer": (0, 'i'),
                 "bwapath": ("bwa", 'i'),
                 "bowtie": ("bowtie", "s"),
                 "check_both": (True, 'b'),
                 "compression": ("bb", 's'),
                 "coorsort": (False, 'b'),
                 "create_log": (False, 'b'),
                 "experiment": ("", 's'),
                 "file_prefix": ("default", 's'),
                 "gatkpath": ("default", 's'),
                 "genome_size": (3.2e9, 'f'),
                 "head": (4, 'i'),
                 "homing": ("CATG", 's'),
                 "intelDeflator": ("default", 's'),
                 "insert_distance": (35, 'i'),
                 "is_slave": (False, 'b'),
                 "keepConsensus": (False, 'b'),
                 "ligation_efficiency": (0.5, 'f'),
                 "logfile": ("default", 's'),
                 "mapped_fraction": (0.83, 'f'),
                 "maxAF": (0.1, 'f'),
                 "MaxPValue": (1e-15, 'f'),
                 "mean_aligned_fraction": (0.9, 'f'),
                 "minBQ": (0, 'i'),
                 "minClustDepth": (10, 'i'),
                 "minCov": (5, 'i'),
                 "minFA": (3, 'i'),
                 "minFM": (3, 'i'),
                 "minFam": (10, 'i'),
                 "minFracAgreed": (0.75, 'f'),
                 "minMQ": (0, 'i'),
                 "minPileupLen": (10, 'i'),
                 "mismatches": (0, 'i'),
                 "mm": (1, 'i'),
                 "num_nucs": (2, 'i'),
                 "on_target": (0.25, 'f'),
                 "opts": ("", 's'),
                 "outBAM": ("default", 's'),
                 "outfile": ("default", 's'),
                 "outfile_handle": (sys.stdout, 'o'),
                 "outTsv": ("default", 's'),
                 "outVCF": ("default", 's'),
                 "overlapLen": (6, 'i'),
                 "parallel": (False, 'b'),
                 "p3Seq": ("default", 's'),
                 "p5Seq": ("default", 's'),
                 "padding": (0, 'i'),
                 "padding_distance": (120, 'i'),
                 "qc_fail": (0.1, 'f'),
                 "readLength": (-1, 'i'),
                 "realigner": ("abra", 's'),
                 "ref": ("default", 's'),
                 "rescue": (False, 'b'),
                 "review_dir": ("default", 's'),
                 "slave_sort_mem": ("768M", 's'),
                 "single_end": (False, 'b'),
                 "sortMem": ("6G", 's'),
                 "threads": (4, 'i'),
                 "uncompressed_bam": (False, 'b')
                 }