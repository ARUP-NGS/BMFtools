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

To run the main program, call the main.py function after installation, or, if installed, run the executable BMFMain.

```bash
python main.py R1.fastq R2.fastq -i BC.fastq -r ${PathToGenomeIndex} --shades --bed ${PathToBedFile}
```

```bash
BMFMain R1.fastq R2.fastq -i BC.fastq -r ${PathToGenomeIndex} --shades --bed ${PathToBedFile}
```

If using a config file, this is greatly simplified. I recommend having multiple config files, one for each analysis type. (e.g., one for FFPE, one for amplicon, one for each bed file, etc.)
A sample configuration file can be found in conf/config.txt. Over time, the number of command-line options has exploded, and a config file is an easy way to keep things consistent.

In that case, one would call thus:

```bash
BMFMain R1.fastq R2.fastq -i BC.fastq --conf ${PathToConfigFile}
```

To use bmftools subcommands, check instructions by executing the following:

```bash
bmftools --help
```

```bash
bmftools <subcommand> --help
```

## Dependencies

Required python packages: Biopython, pysam, pudb, cytoolz, matplotlib, scipy, statsmodels

cutadapt is required for adapter trimming.

numconv is required for conversion to base 64 for PV tags, but that compression is optional.

### Required external tools:
bwa (mem or aln, depending on needs.)

#### Compiler
gcc >= gccv5.0
(Not essential, but gcc5 does offer -floop-unroll-and-jam and a few other neat compiler optimizations.)

#### Adapter Trimming
Cutadapt

#### Indel Realigners
Assembly Based Realigner (abra) (requires bwa)

GATK IndelRealigner


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
FM | Size of family (number of reads sharing barcode.), e.g., "Family Members" | Integer |
FP | Read Passes Filter related to barcoding | For FASTQ: String. Required: "Pass" or "Fail". For BAM: Integer. [0,1] |
ND | Number of Differences in a family of reads from the consensus read. | Integer from Z+ |
NF | ND fraction (mean ND per read in family) | Float |
PV | Phred Values for a read which has saturated the phred scoring system | String, in the form of a comma-joined list of integers. Regex: ASCII|
RA | Realigned due to a failure to map appropriately, either as too small a fraction aligned or a mapping quality of 0. | String: aligned. Current Regex: [a-z]|
RP | Read Pair Position Starts (sorted, separated by a comma) | String. Regex: [GLXYMT0-9.]+:[0-9]+,[GLXYMT0-9.]+[0-9]+ |
SC | Contig Set | String. Regex: [GLXYMT0-9.]+,[GLXYMT0-9]+ |
SF | Soft-Clipped Fraction | Float |
SN | Tags relevant to SNV calling assigned to BAM records. Currently lumped in with SV due to the fact that many are relevant to both.| Comma-separated list of tags. Regex: [A-Z,]+ |
SV | Tags relevant to Structural Variation | Comma-separated list of tags. Regex: [A-Z,]+ |

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


1. Settings Recommendations

    1. The "readPairsPerWrite" parameter can provide great speed improvements.
        1. For my workstation (64GB RAM, 16 threads), the following table indicates that 100 gives me peak performance.
        2. For my cert server (192GB RAM, 24 threads), it looks like 10 might give me peak performance, but more rigorous tests are underway.
    2. For optimal compilation, use the -march flag. BMFTools' setup.py automatically attempts to find that appropriate value for you.

|readPairsPerWrite | time | 
|------|--------------|
| 10 | 867 msec per loop |
| 50 | 851 msec per loop |
| 100 | 850 msec per loop |
| 150 | 853 msec per loop |
|250| 898 msec per loop | 
|500 | 1830 msec per loop |



#TODO (ish):
1. Check and make sure that the bwasw realignment can and does work.
2. Fix the SV tagging - it seems that it broke when I switched to the function call.

## Backlog
1. SNV:
    1. Error Characterization Code
        1. Write database reading and processing.
    3. Info Fields
        2. Add read length to the INFO field (both full read length and read length without Ns)
2. Indels:
    1. Assembly work
    2. Figure out the indel relevance piece.
3. SV:
    1. Finish consensus sequence for intrachromosomal.
    2. Finish writing structural variants to a VCF format
    3. Work on interchromosomal translocations
4. Make sure that the AF filter is being triggered - required 1.0 AF didn't produce expected behavior.
