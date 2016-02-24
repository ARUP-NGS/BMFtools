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


## BMF Tags

These tags are used both in the fastq and the SAM/BAM files.
The only difference between the SAM/BAM tags and the Fastq tags are that the SAM/BAM tags are tab-delimited (as described in sam specifications), while the fastq tags are separated by a delimiter and separated from their values by "=" instead of ":".

Tag | Content | Format |
:----:|:-----|:-----:|
AF | Aligned Fraction aligned (fraction of bases mapped to reference bases, not counting IDSHNP operations. | Float |
DR | Whether the read was sequenced from both strands. Only valid for Loeb-like inline barcodes. | Integer [0, 1] |
FA | Number of reads in Family which Agreed with final sequence at each base | Comma-separated list of integers. Regex: [0-9,]+ |
FM | Size of family (number of reads sharing barcode.), e.g., "Family Members" | Integer |
FP | Read Passes Filter related to barcoding | Integer [0, 1]|
MF | Mate fraction aligned (fraction of bases mapped to reference bases, not counting IDSHNP operations. | Float |
mc | Mate soft-clipped length | Integer |
NC | Number of changed bases in rescued families of reads. | Integer |
NF | Mean number of differences between reads and consensus per read in family | Single-precision floating number |
PV | Phred Values for a read which has saturated the phred scoring system | uint32_t array |
RV | Number of reversed reads in consensus. Only for Loeb-style inline chemistry. | Integer |
SC | Soft-clipped length | Integer |

## Barcoding methods

Essentially, the process is *molecular* demultiplexing.
####Secondary Index Barcoding 
Requires read fastqs and an additional fastq containing barcodes.
> bmftools sdmp
(Secondary-index DeMultiPlex)


####Inline (Loeb-like) barcoding
> bmftools dmp
(DeMultiPlex) 
Barcodes are inline in the start of each read. Because the adapters are enzymatically filled-in, we end with one barcode for each double-stranded template molecule, while the secondary index barcoding ends with 2. This provides better error correction and more accurate diversity quantitation, but the chemistry is much more complicated.

