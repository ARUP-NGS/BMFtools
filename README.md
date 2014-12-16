BMF ( Barcode Manipulation and Factorization)
===================

Suite of tools for barcoded reads taking advantage of PCR redundancy for error reduction/elimination.

## BMF Tags

These tags are used both in the fastq and the SAM/BAM files.
The only difference between the SAM/BAM tags and the Fastq tags are that the SAM/BAM tags are tab-delimited (as described in sam specifications), while the fastq tags are separated by a delimiter and separated from their values by "=" instead of ":".

Tag | Content | Format |
----|-----|-----|
BS | Barcode Sequence | String. Regex: [ATGCN]+ |
FP | Read Passes Filter related to barcoding | String. Required: "Pass" or "Fail" |
FM | Size of family (number of reads sharing barcode.), e.g., "Family Members" | Integer |
BD | Barcode Edit Distance | Integer |

Barcode Edit Distance is 0 for members in a family whose barcode matches the family's exactly. If a rescue step is performed to merge a read with a small number of mismatches due to sequencing errors, this tag will reflect the number of differing characters.

## Barcode Determination methods

####i5/i7 barcoding, nicknamed 'Shades'

Requires read fastqs and an additional fastq containing barcodes.
Faster than using a homing sequence-specified barcode (informatically). More issues with barcode rescues and errors occurring in the auxiliary fastq. Less complicated sample prep.

####Homing sequence-specified regions for barcode.

Using a homing sequence as input for consolidating families of PCR duplicates.
