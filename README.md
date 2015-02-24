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

To use bmftools subcommands, check instructions by executing the following:

```
bmftools --help
```

```
bmftools <subcommand> --help
```

## Dependencies

Required python packages: Biopython, pysam, pudb

### Required external tools:
bwa	aln (for SNV calling)
bwa mem (for SV calling)

### Optional, but recommended:
#### Indel Realigners
Assembly Based Realigner (abra)
GATK IndelRealigner
#### Error Correction Tools
Lighter (Error correction via Bloom Filters)
Reptile (Error correction based on K-mers)


## BMF Tags

These tags are used both in the fastq and the SAM/BAM files.
The only difference between the SAM/BAM tags and the Fastq tags are that the SAM/BAM tags are tab-delimited (as described in sam specifications), while the fastq tags are separated by a delimiter and separated from their values by "=" instead of ":".

Tag | Content | Format |
----|-----|-----|
BS | Barcode Sequence | String. Regex: [ATGCN]+ |
FP | Read Passes Filter related to barcoding | For FASTQ: String. Required: "Pass" or "Fail". For BAM: Integer. [0,1] |
FM | Size of family (number of reads sharing barcode.), e.g., "Family Members" | Integer |
FA | Number of reads in Family which Agreed with final sequence at each base | Comma-separated list of integers. Regex: [0-9,]+ |
BD | Barcode Edit Distance | Integer |
SV | Tags relevant to Structural Variation | Comma-separated list of tags. Regex: [A-Z,]+ |
PV | Phred Values for a read which has saturated the phred scoring system| String, in the form of repr() on a list of integers. Regex: [0-9,\[\]]|
RP | Read Pair Position Starts (sorted, separated by a comma) | String. Regex: [GLXYMT0-9.]+:[0-9]+,[GLXYMT0-9.]+[0-9]+ |
CS | Contig Set | String. Regex: [GLXYMT0-9.]+,[GLXYMT0-9]+ |

## Valid Tags for SV SAM tag

Tag | Meaning |
---- | ----- |
LI | Large Insert - Default cutoff: 1,000,000 min |
MDC | Reads in pair are mapped to different contigs |
MSS | Mapped to Same Strand |
ORB | Only one read in pair mapped to Expected Bed Region |
ORU | One Read Unmapped |
SBI | SBI for having ORB and one of either MDC or LI |
NF | No SV relevance found. |

Barcode Edit Distance is 0 for members in a family whose barcode matches the family's exactly. If a rescue step is performed to merge a read with a small number of mismatches due to sequencing errors, this tag will reflect the number of differing characters.

## Barcode Determination methods

####i5/i7 barcoding, nicknamed 'Shades'

Requires read fastqs and an additional fastq containing barcodes.
Faster than using a homing sequence-specified barcode (informatically). More issues with barcode rescues and errors occurring in the auxiliary fastq. Less complicated sample prep.

####Homing sequence-specified regions for barcode.

Using a homing sequence as input for consolidating families of PCR duplicates.

#Config file

Each line has a set of keys and values. See conf/config.txt for an example.
Most options are available for command-line as well. If an option is set in both a config file and on the command-line, the command-line option clobbers the config setting.

#Changes in BMFTools v0.0.5
1. Removal of standard BMFMain in lieu of the config-based one.
2. Working intrachromosomal translocation detection. (Fast!)
3. Addition of >93 q scores to the read description. This isn't currently used by the variant callers, but it's information which could be used. It does significantly affect the speed of the bmftools dmp step, however.
4. Added filter by bed file to BCVCF. In spite of pysam's supposed ability to pileup over requested reasons, something seems off, so any variants which were called due to pysam's pileup but were outside the bed file are now removed.
5. SNV calling is now in prototypical alpha mode.

#Changes in BMFTools v0.0.5.1
1. Gzipped Fastq's supported.
2. Performance improvements
3. Code now departing from valid python code to cython.

#Changes in BMFTools v.0.0.5.2

1. VCF Info fields for fractions of reads mapped to reverse strand for both alt allele and all reads.
2. VCF Info fields for the mean and standard deviation of base position in read for bases supporting variant call. 
2. Require duplex sequencing an option for variant calling.
3. Adding extra BAM tags for the number of reads in a family supporting the merged read's nucleotide position by position.
4. Moved exceptions to an ErrorHandling file, added PermissionException.
5. Bug fixes for working with gzipped files.


#TODO:
1. Speed up consolidateInexactNumpy (pysam/cStringIO)
2. Finish consensus sequence for intrachromosomal.
3. Finish writing structural variants to a VCF format
4. Work on interchromosomal translocations
5. Error Characterization Code (Start looking at read families differently). Finding a "consensus" sequence for each family, followed by seeing what errors are found at lower family sizes.
7. Consider haplotyping
