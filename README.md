BMF ( Barcode Manipulation and Factorization)
===================

Suite of tools for barcoded reads taking advantage of PCR redundancy for error reduction/elimination.

BMF Tags

These tags are used both in the fastq and the SAM/BAM files.
The only difference between the SAM/BAM tags and the Fastq tags are that the SAM/BAM tags are tab-delimited (as described in sam specifications), while the fastq tags are separated by a delimiter and separated from their values by "=" instead of ":".

Tag | Content | Format |
BS | Barcode Sequence | String. Regex: [ATGCN]+ |
FP | Read Passes Filter related to barcoding | String. Required: "Pass" or "Fail" |
FM | Size of family (number of reads sharing barcode.), e.g., "Family Members" | Integer |

