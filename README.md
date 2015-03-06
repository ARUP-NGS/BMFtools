#BMFTools
###Summary
>BMFTools is a suite of tools for barcoded reads which takes advantage of PCR redundancy for error reduction/elimination.

###What is BMF?
>BMF is a promiscuous acronym with multiple approriate meanings. Primarily:

|By __Reference__| By __Value__ | 
|---------------|-----------|
|_**B**arcoded **M**olecular **F**amilies_ | _**B**arcode **M**anipulation and **F**actorization_ |
######[Cf. Nicklaus Wirth]

===================


## Installation

Run:
```python
python setup.py install
```
You might have an error claiming that README.md is not in dist/. If necessary, copy that file from the distribution base to dist.

## Use

To run the main program, call the main.py function after installation, or, if installed, run the executable BMFMain.

```python
python main.py R1.fastq R2.fastq -i BC.fastq -r ${PathToGenomeIndex} --shades --bed ${PathToBedFile}
```

```python
BMFMain R1.fastq R2.fastq -i BC.fastq -r ${PathToGenomeIndex} --shades --bed ${PathToBedFile}
```

To use bmftools subcommands, check instructions by executing the following:

```python
bmftools --help
```

```python
bmftools <subcommand> --help
```

## Dependencies

Required python packages: Biopython, pysam, pudb

### Required external tools:
bwa (mem or aln, depending on needs.)

### Optional, but recommended:
#### Indel Realigners
Assembly Based Realigner (abra)
GATK IndelRealigner


## BMF Tags

These tags are used both in the fastq and the SAM/BAM files.
The only difference between the SAM/BAM tags and the Fastq tags are that the SAM/BAM tags are tab-delimited (as described in sam specifications), while the fastq tags are separated by a delimiter and separated from their values by "=" instead of ":".

Tag | Content | Format |
----|-----|-----|
BS | Barcode Sequence | String. Regex: [ATGCN]+ |
FP | Read Passes Filter related to barcoding | For FASTQ: String. Required: "Pass" or "Fail". For BAM: Integer. [0,1] |
FM | Size of family (number of reads sharing barcode.), e.g., "Family Members" | Integer |
FA | Number of reads in Family which Agreed with final sequence at each base | Comma-separated list of integers. Regex: [0-9,]+ |
SNV | Tags relevant to SNV calling assigned to BAM records. | Comma-separated list of tags. Regex: [A-Z,]+ |
SV | Tags relevant to Structural Variation | Comma-separated list of tags. Regex: [A-Z,]+ |
PV | Phred Values for a read which has saturated the phred scoring system| String, in the form of repr() on a list of integers in base 85 encoding. Regex: ASCII|
RP | Read Pair Position Starts (sorted, separated by a comma) | String. Regex: [GLXYMT0-9.]+:[0-9]+,[GLXYMT0-9.]+[0-9]+ |
CS | Contig Set | String. Regex: [GLXYMT0-9.]+,[GLXYMT0-9]+ |

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

#Changes in BMFTools v0.0.5:
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

    1. The most import thing is that all of the new features I've been adding finally are debugged and workable.
    1. Performance improvements throughout, but perhaps not entirely relevant.

#Settings Recommendations

    1. The "readPairsPerWrite" parameter can provide great speed improvements. For my workstation, I have what looks like the following:

|readPairsPerWrite | time | 
|------|--------------|
| 10 | 867 msec per loop |
| 50 | 851 msec per loop |
| 100 | 850 msec per loop |
| 150 | 853 msec per loop |
|250| 898 msec per loop | 
|500 | 1830 msec per loop |

Meaning peak performance at 100.


#TODO:
1. SNV:
    1. Consider haplotyping by leveraging reads covering multiple SNPs.
    2. Error Characterization Code (Start looking at read families differently). Finding a "consensus" sequence for each family, followed by seeing what errors are found at lower family sizes.
    3. (Minor) Consider additional SNV tags.
2. Indels:
    1. Work on smaller indels directly in BAM with cigar strings.
    2. Indel realignment might perform better if the "normal" reads are removed, IE, properly-mapped reads without I, D, or S in it.
    3. Try calling freebayes at high ploidy for indels
    4. Try calling scalpel with indel irrelevant reads as normal and relevant reads as abnormal.
    5. Write the DSI SV BAM tag function.
3. SV:
    1. Finish consensus sequence for intrachromosomal.
    2. Write SV tags into a function, call that function during the standard FM/FP/BS tagging.
    3. Finish writing structural variants to a VCF format
    4. Work on interchromosomal translocations
4. Performance:
    1. Continue to implement "map" function for performance gains, especially in BCFastq

5. TODO Backlog/Mostly finished:
    1. Take advantage of PV tags further (done to some extent - moving to backlog)
    2. Instruct pysam to make tags of a specific type.
    3. Note: Something seems wrong with the FA calculation.

1. TODO, lite:
    1. Add control of AlleleAggregateInfo creation (e.g., minFA) to SNVCrawler. **

#Known Issues
1. VCF
    1. NSS INFO field does not work.
    2. MQM, MPF, and MQB have nonsense values.
    3. Need to test NDPS
    4. Strangely-named empty files...???
