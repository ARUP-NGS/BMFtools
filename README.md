#BMFTools
###Summary
>BMFTools is a suite of tools for barcoded reads which takes advantage of PCR redundancy for error reduction/elimination.

###What is BMF?
>BMF is a promiscuous acronym with multiple appropriate meanings. Primarily:

|By __Reference__| By __Value__ |
|---------------|-----------|
|_**B**arcoded **M**olecular **F**amilies_ | _**B**arcode **M**anipulation and **F**actorization_ |
######[Cf. Nicklaus Wirth]

===================


## Installation

Run:
```bash
python setup.py install
```
You might have an error claiming that README.md is not in dist/. If necessary, copy that file from the distribution base to dist.

## Use

To run the main program, call the main.py function after installation, or, if installed, run the executable bmftools main.

```bash
python main.py R1.fastq R2.fastq -i BC.fastq -r ${PathToGenomeIndex} --shades --bed ${PathToBedFile}
```

```bash
bmftools main R1.fastq R2.fastq -i BC.fastq -r ${PathToGenomeIndex} --shades --bed ${PathToBedFile}
```

If using a config file, this is greatly simplified. I recommend having multiple config files, one for each analysis type. (e.g., one for FFPE, one for amplicon, one for each bed file, etc.)
A sample configuration file can be found in conf/config.txt. Over time, the number of command-line options has exploded, and a config file is an easy way to keep things consistent.

In that case, one would call thus:

```bash
bmftools main R1.fastq R2.fastq -i BC.fastq --conf ${PathToConfigFile}
```

To use bmftools subcommands, check instructions by executing the following:

```bash
bmftools --help
```

```bash
bmftools <subcommand> --help
```

## Dependencies

Required python packages: Biopython, pysam, pudb, cytoolz, matplotlib, scipy, statsmodels, numconv

cutadapt is required for adapter trimming.

numconv is required for conversion to base 64 for PV tags, but that compression is optional.

re2 (Google) is a fast alternative to standard python re. BMFTools will attempt to load the re2 as re. Failing that, it will fall back to the standard library.

Some of vcflib's tools are used, although vcflib (key-word argument for these calls) can be set to False to do it manually on the shell.

### Required external tools:

#### Utilities set
samtools >= 1.1
bamleftalign (from FreeBayes)

#### Aligners

bwa >= 0.7.10 (mem, aln, bwasw)

bowtie/bowtie2


#### Compiler
gcc >= gccv5.0
(Not essential, but gcc5 does offer -floop-unroll-and-jam and a few other neat compiler optimizations.)

#### Adapter Trimming
Cutadapt >= 1.7

#### Indel Realigners
ABRA >= 0.85 (Assembly Based Realigner, which in turn requires bwa)

GATK >= 1.6, for its IndelRealigner

#### Depth Calculations
[FastDepthOfCoverage](https://github.com/ARUP-NGS/Pipeline/blob/master/src/main/java/util/coverage/CovCalcApp.java) [Now included in java/ folder as FastDOC.java]
Alternatively, you may use GATK.


## BMF Tags

These tags are used both in the fastq and the SAM/BAM files.
The only difference between the SAM/BAM tags and the Fastq tags are that the SAM/BAM tags are tab-delimited (as described in sam specifications), while the fastq tags are separated by a delimiter and separated from their values by "=" instead of ":".

Tag | Content | Format |
----|-----|-----|
AF | Aligned Fraction | Float |
BS | Barcode Sequence | String. Regex: [ATGCN]+ |
CC | Cluster Count | Integer |
DP | Discordant Read Pair information. | String |
FA | Number of reads in Family which Agreed with final sequence at each base | Comma-separated list of integers. Regex: [0-9,]+ |
FF | Minimum value for Fraction of FA for a base in a condensed read. | Float [0-1] | Regex: [0-9\.]+ |
PF | Ratio of minimum and maximum PV values for bases in a condensed read. | Float [0-1] | Regex: [0-9\.]+ |
FM | Size of family (number of reads sharing barcode.), e.g., "Family Members" | Integer |
FP | Read Passes Filter related to barcoding | For FASTQ: String. Required: "Pass" or "Fail". For BAM: Integer. [0,1] |
ND | Number of Differences in a family of reads from the consensus read. | Integer from Z+ |
NF | ND fraction (mean ND per read in family) | Float |
PV | Phred Values for a read which has saturated the phred scoring system | String, in the form of a comma-joined list of integers. Regex: ASCII|
RA | Realigned due to a failure to map appropriately, either as too small a fraction aligned or a mapping quality of 0. | String: aligned. Current Regex: [a-z]|
RP | Read Pair Position Starts (sorted, separated by a comma) | String. Regex: [GLXYMT0-9\.]+:[0-9]+,[GLXYMT0-9\.]+[0-9]+ |
SC | Contig Set | String. Regex: [GLXYMT0-9\.]+,[GLXYMT0-9]+ |
SF | Soft-Clipped Fraction | Float |
SV | Tags relevant to Structural Variation | Comma-separated list of tags. Regex: [A-Z,]+ |

## Read Pair Merging Tags

These are only used for merging read pairs.

Tag | Content | Format |
----|-----|-----|
DG | Discordant positions in merged pair, genomic coordinates. | String. Regex: [0-9,]+ |
DR | Discordant read positions in merged pair, genomic coordinates. | String. Regex: [0-9,]+ |
MA | Indices for read positions which agreed during merging. | String. Regex: [0-9,]+ |
mp | Original Mate Position | Integer |
om | Original Mapping Quality | Integer |
op | Original Position | Integer |
ot | Original Template Langth | Integer |
PM | Indices for read positions which have been merged | String. Regex: [0-9,]+ |

## Valid Tags for SV SAM tag

######Note: SNV tags have been lumped in with SV tags for the time being.

Tag | Meaning |
---- | ----- |
LI | Large Insert - Default cutoff: 1,000,000 min |
MDC | Reads in pair are mapped to different contigs |
MSS | Mapped to Same Strand |
ORB | Only one read in pair mapped to Expected Bed Region |
ORU | One Read Unmapped |
ORS | One Read Is Soft-Clipped |
DRP | Duplex Read Pair |
DSI | Duplex Supported Insertion |
DSD | Duplex Supported Deletion |
NF | No SV relevance found. |


## Barcode Determination methods

####i5/i7 barcoding, nicknamed 'Shades'

Requires read fastqs and an additional fastq containing barcodes.
Faster than using a homing sequence-specified barcode (informatically). More issues with barcode rescues and errors occurring in the auxiliary fastq. Less complicated sample prep.

####Homing sequence-specified regions for barcode.

Using a homing sequence as input for consolidating families of PCR duplicates.

#Config file

Each line has a set of keys and values. See conf/config.txt for an example.
Most options are available for command-line as well. If an option is set in both a config file and on the command-line, the command-line option clobbers the config setting.

1. Changes in BMFTools v0.0.5:
    1. Removal of standard BMFMain in lieu of the config-based one.
    2. Working intrachromosomal translocation detection. (Fast!)
    3. Addition of >93 q scores to the read description. This isn't currently used by the variant callers, but it's information which could be used. It does significantly affect the speed of the bmftools dmp step, however.
    4. Added filter by bed file to BCVCF. In spite of pysam's supposed ability to pileup over requested reasons, something seems off, so any variants which were called due to pysam's pileup but were outside the bed file are now removed.
    5. SNV calling is now in prototypical alpha mode.

1. Changes in BMFTools v0.0.5.1:
    1. Gzipped Fastq's supported.
    2. Performance improvements
    3. Code now departing from valid python code to cython.

2. Changes in BMFTools v0.0.5.2:

    1. VCF Info fields for fractions of reads mapped to reverse strand for both alt allele and all reads.
    2. VCF Info fields for the mean and standard deviation of base position in read for bases supporting variant call.
    2. Require duplex sequencing an option for variant calling.
    3. Adding extra BAM tags for the number of reads in a family supporting the merged read's nucleotide position by position (FA tag).
    4. Moved exceptions to an ErrorHandling file, added PermissionException.
    5. Bug fixes for working with gzipped files.

3. Changes in BMFTools v0.0.5.3:

    1. Sped-up, cythonized, and numpied version of fastq consolidation fully functional and fast.
    2. Compiler directives.
    3. More type definitions.
    4. Optionally encoding the summed phred scores from demultiplexed families in base 85 format to save space.
    5. Removed the "slow" form of fastq consolidation.

4. Changes in BMFTools v0.0.6.0:

    1. The most important thing is that all of the new features I've been adding finally are debugged and workable.
    1. Performance improvements throughout, but perhaps not entirely relevant.

5. Changes in BMFTools v0.0.6.1:

    1. Created a pFastqProxy object to cause pysam's FastqProxy object's information to persist.
    2. Compiler optimizations

6. Changes in BMFTools v0.0.7.0

    1. Created a pPileupColumn object to cause pysam's PileupColumn object's information to persist. (Sound familiar?)
    2. All quality scores of "2" are replaced by "0" in demultiplexing because they mean nothing.
    3. Replaced a vectorized function calling a dictionary into a list comprehension of that dictionary. (It's faster)
    4. Faster string operations in BCFastq

7. Changes in BMFTools v0.0.7.1

    1. Fixed NSS field.
    2. Verified that NDPS field is working.
    3. MQM and MQB are now working.
    4. PV tags used for all reads now, making compatibility a little easier.

8. Changes in BMFTools v0.0.7.2

    1. Probabilistic quantitation of AAF given observations.
    2. Optional filter for FFPE data for removing deamination frequencies due to formalin fixation.
    3. Added amplicon filtering for mispriming, both in variant caller and in a pre-processing step.
    4. Added a test for allelic imbalance.
    5. Removed "N"s from variant calls.

9. Changes in BMFTools v0.0.7.3

    1. Removing discordant read pairs from pileups for SNV calls.
    2. Added sort memory options.
    3. minAF filter (Aligned Fraction)
    4. Added bwasw realignment for reads failing a given filter.
    5. Added QC measurements and steps.
    6. Additional info fields, BAM tags
    7. Overrode __str__ for objects which had a ToString function.
    8. Performance enhancements.
    9. Further work on indel calls.

10. Changes in BMFTools v0.0.7.4

    1. Fixed NDP calculation.
    2. Fixed AF calculation.
    3. Heavily cythonized
    4. Memoization
    5. Fixed some indel work
    6. Fixed some VCF comparison issues.

11. Changes in BMFTools v0.0.7.5
    1. Parallelized variant-calling
    2. General Popen Dispatcher framework.
    3. Cursory work on tumor/normal pairs for SNVs.
    4. Further indel work.
    5. Optional bwasw realignment for all reads with AF <= minAF

12. Changes in BMFTools v.0.1.0.0beta and v0.1.0.1beta
    0. Re-coded a lot of the under-the-hood stuff.
        1. Fastq consolidation.
        2. Fastq marking
        3. Bam tagging
        4. Alignment
        5. Bam record processing.
        6. Addition of the header line with RG:default in one pass (no call to Picard)
        7. Re-wrote the SV tagging.
        8. Completely re-did barcode rescue.
    1. Wrote the pairedBarcodeTagging in C++ with bamtools API. Defaults to running in cython after piping.
    1. Support for more single-end analysis.
    2. Unit tests.
    3. Addition of optional "head" parameter for salting molecular barcodes.
    4. Indel work.
        1. Toy indel caller that requires that both reads in a pair support it and that the Shannon entropy be above a threshold. That being said, it's not that great, and we still have to try to call variants in low complexity regions.
    5. Uniqueness/mappability work.
        1. Read binning based on an expanded hashmap for O(1) inexact string comparison.
        2. It would be fast if all of our reads were on-target...
    6. A number of optimizations, including some string comparisons as integers.
    8. Removing old/dead code.

13. Changesin BMFTools v0.1.0.2beta
    1. Complete restructure of main program, argument parsing, and error handling.
    2. Moved to a global config/started review directory infrastructure.
    3. Faster failure of pileups where all reads failed filters.
    4. String handling/other under-the-hood.
    5. Chapman now controls all global arguments.


## Style guidelines

### Modified PEP8
Since cython is so different from python, some of the PEP8 rules are (inappropriately) being triggered by this code when pep8 is run on it.
We have decided to follow a modified PEP8, where all PEP8 complaints are considered valid except for "whitespace around operator" for cdefs as appropriate (such as ndarray[float, ndim=1]) and "module level import" not at top of file.
The second is because pep8 doesn't realize that cimports are still imports.
