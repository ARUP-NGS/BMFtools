BMF ( Barcode Manipulation and Factorization)
===================

Suite of tools for barcoded reads taking advantage of PCR redundancy for error reduction/elimination.

## Installation

Run:
```
python setup.py install
```
You might have an error claiming that README.md is not in dist/. If necessary, copy that file from the distribution base to dist.

## Use

To run the main program, call the main.py function after installation, or, if installed, run the executable BMFMain.

```
python main.py R1.fastq R2.fastq -i BC.fastq -r ${PathToGenomeIndex} --shades --bed ${PathToBedFile}
```

```
BMFMain R1.fastq R2.fastq -i BC.fastq -r ${PathToGenomeIndex} --shades --bed ${PathToBedFile}
```

## BMF Tags

These tags are used both in the fastq and the SAM/BAM files.
The only difference between the SAM/BAM tags and the Fastq tags are that the SAM/BAM tags are tab-delimited (as described in sam specifications), while the fastq tags are separated by a delimiter and separated from their values by "=" instead of ":".

Tag | Content | Format |
----|-----|-----|
BS | Barcode Sequence | String. Regex: [ATGCN]+ |
FP | Read Passes Filter related to barcoding | String. Required: "Pass" or "Fail" |
FM | Size of family (number of reads sharing barcode.), e.g., "Family Members" | Integer |
BD | Barcode Edit Distance | Integer |
SV | Tags relevant to Structural Variation | Comma-separated list of tags. Regex: [A-Z,]+ |
PV | Phred Values for a read which has saturated the phred scoring system| String, in the form of repr() on a list of integers. Regex: [0-9,\[\] ]| 

## Valid Tags for SV SAM tag

Tag | Meaning |
---- | ----- |
LI | Large Insert - Default cutoff: 1,000,000 min |
MDC | Reads in pair are mapped to different contigs |
MSS | Mapped to Same Strand |
ORB | Only one read in pair mapped to Expected Bed Region |
ORU | One Read Unmapped |
SBI | SBI for having ORB and one of either MDC or LI |
NF | No SV relevance found found. |

Barcode Edit Distance is 0 for members in a family whose barcode matches the family's exactly. If a rescue step is performed to merge a read with a small number of mismatches due to sequencing errors, this tag will reflect the number of differing characters.

## Barcode Determination methods

####i5/i7 barcoding, nicknamed 'Shades'

Requires read fastqs and an additional fastq containing barcodes.
Faster than using a homing sequence-specified barcode (informatically). More issues with barcode rescues and errors occurring in the auxiliary fastq. Less complicated sample prep.

####Homing sequence-specified regions for barcode.

Using a homing sequence as input for consolidating families of PCR duplicates.

# Config JSON

The run JSON has two main fields: "Config", which contains the system variables and general preferences, and "Analysis", which controls the the run protocol. 
Previously, there were two config files, but they have been consolidated.
