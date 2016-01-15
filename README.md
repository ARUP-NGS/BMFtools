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
DR | Whether the read was sequenced from both strands. Only valid for Loeb-like inline barcodes. | Integer [0, 1] |
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
