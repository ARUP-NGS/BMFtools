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

```bash
git clone https://github.com/ARUP-NGS/BMFtools --recursive
cd BMFtools
make
```
## Use


```bash
bmftools
```

```bash
bmftools <--help/-h>
```

```bash
bmftools <subcommand> <-h>
```

## Dependencies

cutadapt is required for adapter trimming.

pigz is required for compression/decompression

### Required external tools:

#### Utilities set
samtools >= 1.2

#### Aligners

bwa >= 0.7.10 (mem, aln, bwasw)

#### Adapter Trimming
Cutadapt >= 1.7

## BMF Tags

These tags are used both in the fastq and the SAM/BAM files.
The only difference between the SAM/BAM tags and the Fastq tags are that the SAM/BAM tags are tab-delimited (as described in sam specifications), while the fastq tags are separated by a delimiter and separated from their values by "=" instead of ":".

Tag | Content | Format |
:----:|:-----|:-----:|
FA | Number of reads in Family which Agreed with final sequence at each base | Comma-separated list of integers. Regex: [0-9,]+ |
FM | Size of family (number of reads sharing barcode.), e.g., "Family Members" | Integer |
FP | Read Passes Filter related to barcoding | Integer [0, 1]|
NC | Number of changed bases in rescued families of reads. | Integer |
NF | Mean number of differences between reads and consensus per read in family | Single-precision floating number |
PV | Phred Values for a read which has saturated the phred scoring system | uint32_t array|
RV | Number of reversed reads in consensus. Only for Loeb-style inline chemistry. | Integer |

## Read Pair Merging Tags

These are only used for merging read pairs.

Tag | Content | Format |
:----:|:-----|:-----:|
DG | Discordant positions in merged pair, genomic coordinates. | String. Regex: [0-9,]+ |
DR | Discordant read positions in merged pair, genomic coordinates. | String. Regex: [0-9,]+ |
MA | Indices for read positions which agreed during merging. | String. Regex: [0-9,]+ |
mp | Original Mate Position | Integer |
om | Original Mapping Quality | Integer |
op | Original Position | Integer |
ot | Original Template Langth | Integer |
PM | Indices for read positions which have been merged | String. Regex: [0-9,]+ |

## Barcoding methods

####i5/i7 barcoding, nicknamed 'Shades'

Requires read fastqs and an additional fastq containing barcodes.
Faster than using a homing sequence-specified barcode (informatically). More issues with barcode rescues and errors occurring in the auxiliary fastq. Less complicated sample prep.

####Inline (Loeb-like) barcoding

Barcodes are inline in the start of each read. This information is removed, added to a fastq comment, and then used in final hashmap-powered molecular demultiplexing.

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

13. Changes in BMFTools v0.1.0.2beta
    1. Complete restructure of main program, argument parsing, and error handling.
    2. Moved to a global config/started review directory infrastructure.
    3. Faster failure of pileups where all reads failed filters.
    4. String handling/other under-the-hood.
    5. Chapman now controls all global arguments.

14. Changes in BMFtools v0.2.0
    1. Rewritten in C.
    2. Positional rescue for barcoded mismatches implemented.
    3. Error rate calculation utilities.
    4. Fork on samtools sort for rescue.

