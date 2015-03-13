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

cutadapt is required for adapter trimming.

numconv is required for conversion to base 64 for PV tags, but that compression is optional.

### Required external tools:
bwa (mem or aln, depending on needs.)

#### Compiler
gcc >= gccv4.5

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
BS | Barcode Sequence | String. Regex: [ATGCN]+ |
CS | Contig Set | String. Regex: [GLXYMT0-9.]+,[GLXYMT0-9]+ |
FA | Number of reads in Family which Agreed with final sequence at each base | Comma-separated list of integers. Regex: [0-9,]+ |
FM | Size of family (number of reads sharing barcode.), e.g., "Family Members" | Integer |
FP | Read Passes Filter related to barcoding | For FASTQ: String. Required: "Pass" or "Fail". For BAM: Integer. [0,1] |
ND | Number of Differences in a family of reads from the consensus read. | Integer from Z+ |
NF | ND fraction (mean ND per read in family) | Float |
PV | Phred Values for a read which has saturated the phred scoring system| String, in the form of repr() on a list of integers in base 85 encoding. Regex: ASCII|
RP | Read Pair Position Starts (sorted, separated by a comma) | String. Regex: [GLXYMT0-9.]+:[0-9]+,[GLXYMT0-9.]+[0-9]+ |
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

    1. The most import thing is that all of the new features I've been adding finally are debugged and workable.
    1. Performance improvements throughout, but perhaps not entirely relevant.

5. Changes in BMFTools v0.0.6.1:

    1. Created a pFastqProxy object to cause pysam's FastqProxy object's information to persist.
    2. Compiler optimizations

6. Changes in BMFTools v0.0.7.0

    1. Created a pPileupColumn object to cause pysam's PileupColumn object's information to persist. (Sound familiar?)
    2. All quality scores of "2" are replaced by "0" in demultiplexing because they mean nothing.
    3. Replaced a vectorized function calling a dictionary into a list comprehension of that dictionary. (It's faster)
    4. Faster string operations in BCFastq

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
1. SNV:
    0. Just fix everything (See Known Issues)
    1. Add INFO fields for the new NF tag (mean, max, SD)
    1. Consider haplotyping by leveraging reads covering multiple SNPs.
    2. Error Characterization Code (Start looking at read families differently). Finding a "consensus" sequence for each family, followed by seeing what errors are found at lower family sizes.
    3. Add tile/x pos/y pos information to error characterization code, along with cluster count. (Parse from the Stats xml from bcl2fastq2)
    4. Remove BS INFO field, change how the FILTER requiring both strands.
2. Indels:
    1. Work on smaller indels directly in BAM with cigar strings. (Some progress on it, but it has to wait until SNVs are fully-solved)
    2. Indel realignment might perform better if the "normal" reads are removed, IE, properly-mapped reads without I, D, or S in it.
    3. FreeBayes with a longer --haplotype-length, demultiplexing first, and a high ploidy + pre-filtering should get us what we want.
    4. And, perhaps we need something like Scalpel for larger indels.
3. SV:
    1. Finish consensus sequence for intrachromosomal.
    2. Finish writing structural variants to a VCF format
    3. Work on interchromosomal translocations
    4. K-mer based SV approaches?

1. TODO, lite:
    1. Add control of AlleleAggregateInfo creation (e.g., minFA) to SNVCrawler. **

#Known Issues
1. VCF
    1. NSS INFO field does not seem to work.
    2. MQM, MPF, and MQB have nonsense values.
    3. Need to test NDPS

