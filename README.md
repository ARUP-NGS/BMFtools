#BMFTools
###Summary
>BMFtools (_**B**arcoded **M**olecular **F**amilies_tools) is a suite of tools for barcoded reads which takes advantage of PCR redundancy for error reduction/elimination.

===================


### Installation

```bash
git clone https://github.com/ARUP-NGS/BMFtools --recursive
cd BMFtools
make
```

### Tools

Name | Use |
:---:|:----|
bmftools cap| Postprocess a tagged BAM for BMF-agnostic tools.|
bmftools depth| Calculates depth of coverage over a set of bed intervals.|
bmftools dmp| Demultiplex inline barcoded experiments.|
bmftools err| Calculate error rates based on cycle, base call, and quality score.|
bmftools famstats| Calculate family size statistics for a bam alignment file.|
bmftools filter| Filter or split a bam file by a set of filters.|
bmftools hashdmp| Demultiplex inline barcoded experiments that have already been marked.|
bmftools mark| Add unclipped start position as annotation for both read and mate.|
bmftools rsq| Rescue bmf-sorted or ucs-sorted bam alignments.|
bmftools sdmp| Demultiplex secondary-index barcoded experiments.|
bmftools sort| Sort for bam rescue.|
bmftools stack| A maximally-permissive yet statistically-thorough variant caller using molecular barcode metadata.|
bmftools target| Calculates on-target rate.|
bmftools vet| Curate variant calls from another variant caller (.vcf) and a bam alignment.|


### Use

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

Essentially, the process is *molecular* **d**e**m**ulti**p**lexing.

####Secondary Index Barcoding 
Requires read fastqs and an additional fastq containing barcodes.
> bmftools sdmp

(Secondary-index DeMultiPlex)

Barcodes are contained in the additional (typically "Read 2" of 3) fastq, with optional "salting" from the tthe start of the insert reads.
Note: It is **STRONGLY** recommended that for the secondary-index chemistry that you mask adapter sequence in the molecular barcode reads.
This will N those bases, failing the reads, If not, you'll get a ton of garbage!

####Inline (Loeb-like) barcoding
Requires only read fastqs, as the barcodes are inline.
> bmftools dmp

(DeMultiPlex)

Barcodes are inline in the start of each read. Because the adapters are enzymatically filled-in, we end with one barcode for each double-stranded template molecule, while the secondary index barcoding ends with 2. This provides better error correction and more accurate diversity quantitation, but the chemistry is more complicated.

