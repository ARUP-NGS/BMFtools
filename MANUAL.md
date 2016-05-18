#BMFtools
##Summary
BMFtools (**B**arcoded **M**olecular **F**amilies tools) is a suite of tools for error correction and precise quantitation using reads with molecular barcodes.
Reads which are PCR duplicates of original template molecules are **molecularly** demultiplexed into single unique observations for each sequenced founding template molecule.
Accessory tools provide postprocessing, filtering, quality control, and summary statistics.

All command-line options are available from the command-line as follows:

For a list of subcommands:
```bash
bmftools <--help/-h>
```

For usage for a subcommand:
```bash
bmftools <subcommand> <-h>
```

##Synopsis
<p>

`bmftools cap -f 0.8 -c 200 input.bam output.bam`

<p>

`bmftools depth -H coverage_uniformity.hist.txt -b capture.bed -Q1 input.bam > coverage.bed.txt`

<p>

`bmftools dmp -p 12 -d -f output_prefix -l <barcode_length> input_R1.fastq.gz input_r2.fastq.gz`

<p>

`bmftools err fm reference.fasta input.srt.bam > err_by_fm.txt`

<p>

`bmftools err main -c cycle_err.txt -n base_call_err.txt -g global_err.txt -o rescaled_qualities.txt reference.fasta input.srt.bam`

<p>

`bmftools err region -b regions.bed -o err_by_region.txt reference.fasta input.srt.bam`

<p>

`bmftools famstats fm input.bam > famstats.txt`

<p>

`bmftools famstats frac 2 input.bam > frac_fm_ge2.txt`

<p>

`bmftools filter -s2 input.bam output.bam`

<p>

`bmftools mark input.bam output.bam`

<p>

`bmftools rsq -u -f tmp.fastq input.bam tmp.out.bam`

<p>

`bmftools sdmp -p 12 -s 1 -o tmp_prefix -d -f output_prefix -i index.fastq.gz read1.fastq.gz read2.fastq.gz`

<p>

`bmftools sort -T tmp_prefix -k ucs input.bam output.bam`

<p>

`bmftools stack --ref reference.fasta -o output.vcf -b capture.bed --min-family-size 3 tumor.bam normal.bam`

<p>

`bmftools target -b capture.bed input.bam`

<p>

`bmftools vet -o output.vcf -b capture.bed --min-family-size 3 input.bcf input.bam`

##Commands and options

### Core Functionality

####<b>dmp</b>
  Description:
  > Performs molecular demutiplexing of inline barcoded fastq data.
  > The homing sequence is a sequence of bases marking the end of the random nucleotides
  > which make up the barcode. This is required to ensure that chemistr has worked as expected
  > and to identify the barcode length, which is used for additional entropy.

  > To achieve linear performance with arbitrarily large datasets, an initial marking step subsets the reads by the first
  > few nucleotides in the barcode. The more of these are used, the lower the RAM requirements but the more temporary files are written.
  > This is controlled by the -n option.

  Usage: `bmftools dmp <options> input_R1.fastq.gz input_R2.fastq.gz`

  Options:

    > -D:    Skip final consolidation and only create temporary marked subset files.
    > -n:    Number of bases in the beginning of the barcode to use for subsetting.
    > -S:    Run in single-end mode. (ignores read 2)
    > -=:    Emit output to stdout, interleaved if paired-end, instead of writing to disk.
    > -s:    Homing sequence. REQUIRED.
    > -l:    Barcode length. REQUIRED. For variable-length barcodes, this signifies the minimum length.
    > -v:    Maximum barcode length. Only needed for variable-length barcodes.
    > -o:    Temporary file basename. Defaults to a random string variation on the input filename.
    > -t:    Reads with a homopolymer of threshold <parameter> length or greater are marked as QC fail. Default: 10.
    > -m:    Skip first <parameter> bases at the beginning of each read for use in barcode due to their high error rates.
    > -p:    Number of threads to use for dmp step.
    > -f:    Sets final fastq prefix. Final filenames will be <parameter>.R[12].fq if uncompressed, <parameter>.R[12].fq.gz if compressed. Ignored if -= is set.
    > -r:    Path to text file with rescaled quality scores. Used for rescaling quality scores during dmp. Only used if provided.
    > -z:    Flag to write gzip-compressed output.
    > -T:    Write temporary fastq files with gzip compression level <parameter>. Defaults to transparent gzip files (zlib >= 1.2.5) or uncompressed (zlib < 1.2.5).
    > -g:    Gzip compression parameter when writing gzip-compressed output. Default: 1.
    > -u:    Notification interval. Log each <parameter> sets of reads processed during the initial marking step. Default: 1000000.
    > -w:    Leave temporary files.
    > -h/-?: Print help menu.


####<b>sdmp</b>
  Description:
  > Performs molecular demutiplexing of secondary index barcoded fastq data.

  Usage: `bmftools sdmp <options> input_R1.fastq.gz input_R2.fastq.gz`

  Options:

    > -i:    Path to index fastq. REQUIRED.
    > -n:    Number of bases in the beginning of the barcode to use for subsetting.
    > -D:    Skip final consolidation and only create temporary marked subset files.
    > -S:    Run in single-end mode. (ignores read 2)
    > -s:    Number of bases from the beginning of each read to use to "salt" the barcode for additional entropy.
    > -o:    Temporary file basename. Defaults to a random string variation on the input filename.
    > -t:    Reads with a homopolymer of threshold <parameter> length or greater are marked as QC fail. Default: 10.
    > -m:    Skip first <parameter> bases at the beginning of each read for use in barcode salting due to their high error rates.
    > -p:    Number of threads to use for dmp step.
    > -=:    Emit output to stdout, interleaved if paired-end, instead of writing to disk.
    > -f:    Sets final fastq prefix. Final filenames will be <parameter>.R[12].fq if uncompressed, <parameter>.R[12].fq.gz if compressed. Ignored if -= is set.
    > -r:    Path to text file with rescaled quality scores. Used for rescaling quality scores during dmp. Only used if provided.
    > -z:    Flag to write gzip-compressed output.
    > -T:    Write temporary fastq files with gzip compression level <parameter>. Defaults to transparent gzip files (zlib >= 1.2.5) or uncompressed (zlib < 1.2.5).
    > -g:    Gzip compression parameter when writing gzip-compressed output. Default: 1.
    > -u:    Notification interval. Log each <parameter> sets of reads processed during the initial marking step. Default: 1000000.
    > -w:    Leave temporary files.
    > -h/-?: Print help menu.


####<b>rsq</b>
  Description:
  > Uses positional information to rescue reads into proper families in cases were there were errors in the barcode.
  > In preprocessing, bam file is sorted to group reads sharing an "alignment signature" together, where an alignment signature consists of position, orientation, and read length for a read (and its mate, for paired-end data).
  > In a process analogous to samtools rmdup, reads sharing an alignment signature with a barcode hamming distance
  > below a given threshold are considered to have originated from the same template molecule and are consolidted into one observation.
  > For paired-end data, requires pre-processing by bmftools mark to add mate information as auxiliary tags.
  > For all data, requires sorting by alignment signature (see bmftools sort).
  > Because rescued read families may have changed base calls, these reads written to a temporary fastq
  > and realigned. In addition, this regenerates all of the secondary and supplementary alignments.

  Usage: `bmftools rsq <options> input_R1.srt.bam output.bam`

  Options:

    > -f:    Path to temporary fastq
    > -S:    Flag to perform rescue for single-end experiments.
    > -u:    Flag to use unclipped start rather than position in alignment signature. Requires `-k ucs` in bmftools sort.
    > -s:    Flag to write reads with supplementary alignments to the temporary fastq to regenerate secondary/supplementary reads after rescue.
    > -l:    Output bam compression level.
    > -t:    Mismatch limit. Default: 2.
    > -h/-?: Print help menu.

### Analysis

####<b>stack</b>
    Description:
    > A maximally-permissive single nucleotide variant caller for matched sample pairs using molecular barcode metadata analogous to samtools mpileup.
    > Reads with the same read name are merged into a single observation for each pileup, with p-values merged according to Fisher's method.
    > Passing calls are marked BMF_PASS.
    > Passing calls in the tumor but not the normal are marked as SOMATIC.

    Options:
    > -R, --refpath:               Path to fasta reference. REQUIRED.
    > -b, --bed-path:              Path to bed file for anaylsis. REQUIRED.
    > -p, --padding:               Number of bases around each region to pad in calling variants.

    > -o, --outpath:               Path to output [bv]cf. Defaults to stdout.
    > -c, --min-count:             Minimum number of passing observations to pass a variant.
    > -s, --min-family-size:       Minimum family size required to pass an observation.
    > -f, --min-fraction-agreed:   Minimum fraction of family members agreed on a base call.
    > -v, --min-phred-quality:     Minimum PV tag value required.
    > -a, --min-family-agreed:     Minimum number of reads in a family agreed on a base call.
    > -m, --min-mapping-quality:   Minimum mapping quality required for inclusion of a read.

    > -2, --skip-secondary:        Skip secondary alignments.
    > -S, --skip-supplementary:    Skip supplementary alignments.
    > -q, --skip-qc-fail:          Skip reads marked as QC fail.
    > -r, --skip-duplicates:       Skip reads marked as duplicates.
    > -B, --emit-bcf-format:       Emit bcf-formatted output instead of vcf.

    TODO: Fill in details on these tags.
    VCF Header Fields:
    1. BMF_PASS
    2. ADP
    3. ADPO
    3. ADPD
    4. ADPR
    5. RVF
    6. QSS
    7. AMBIG
    8. SOMATIC_PV
    9. SOMATIC_CALL
    10. SOMATIC
    11. FR_FAILED
    12. FM_FAILED
    13. FP_FAILED
    14. AF_FAILED
    15. MQ_FAILED
    16. IMPROPER
    17. OVERLAP

####<b>vet</b>
    Description:
    > Curates variant calls from a bcf file and an associated, indexed bam file.

    Options:
    > -b, --bed-path:              Path to bed file for anaylsis. REQUIRED.
    > -p, --padding:               Number of bases around each region to pad in calling variants.

    > -o, --outpath:               Path to output [bv]cf. Defaults to stdout.
    > -c, --min-count:             Minimum number of passing observations to pass a variant.
    > -s, --min-family-size:       Minimum family size required to pass an observation.
    > -f, --min-fraction-agreed:   Minimum fraction of family members agreed on a base call.
    > -v, --min-phred-quality:     Minimum PV tag value required.
    > -a, --min-family-agreed:     Minimum number of reads in a family agreed on a base call.
    > -m, --min-mapping-quality:   Minimum mapping quality required for inclusion of a read.

    > -2, --skip-secondary:        Skip secondary alignments.
    > -S, --skip-supplementary:    Skip supplementary alignments.
    > -q, --skip-qc-fail:          Skip reads marked as QC fail.
    > -F, --skip-recommended:      Skip secondary, supplementary, and PCR duplicates.
    > -B, --emit-bcf-format:       Emit bcf-formatted output instead of vcf.
    > -w, --write-outside-bed:     Write variants outside of bed region unmodified rather than removing.

    TODO: Fill in details on these tags.
    VCF Header Fields:
    1. BMF_VET
    2. BMF_UNIOBS
    3. BMF_DUPLEX
    4. BMF_FAIL
    5. DUPLEX_DEPTH
    6. DISC_OVERLAP
    7. OVERLAP


### Manipulation

####<b>cap</b>
  Description:
  > Caps quality scores using barcode metadata to facilitate working with barcode-agnostic tools.

  Usage: `bmftools cap <options> input_R1.srt.bam output.bam`

  Options:

    > -l:    Set output compression level. Default: 6.
    > -t:    Set phred score to which to set passing base qualities. Default: 93 ('~').
    > -m:    Set minFM required to pass reads. Default: 0.
    > -f:    Minimum fraction of reads in a family supporting a base call for inclusion. Default: 1.0.
    > -c:    Set minimum calculated phred score to not mask a base call. Default: 0. 
    > -d:    Flag to only mask failing base scores as '#'/2, not modifying passing quality scores.
    > -h/-?: Print help menu.

####<b>filter</b>
  Description:
  > Filters or splits a bam file. In filter mode, only passing reads are output. In split mode,
  > emits passing reads to one file and failing reads to another.

  Usage: `bmftools filter <options> input_R1.srt.bam output.bam`

  Options:

    > -l:    Set output compression level. Default: 6.
    > -a:    Read pairs without one read with an aligned fraction above <parameter> are failed.
    > -m:    Minimum mapping quality.
    > -F:    Fail all reads with any bits in <parameter> set.
    > -f:    Fail all reads without all bits in <parameter> set.
    > -b:    Require reads be within the region defined by the bed file at <parameter>.
    > -P:    Number of bases around the bed file with which to pad.
    > -r:    If set, write failing reads to bam at <parameter>.
    > -v:    Invert pass/fail. (Analogous to grep.)

####<b>depth</b>
  Description:
  > Creates a bed file of coverage depths for both raw and collapsed read families over a capture region of interest.
  > Requires an indexed bam.

  Usage: `bmftools depth <options> input_R1.srt.bam`

  Options:

    > -o:    Write coverage bed to <path> instead of stdout.
    > -H:    Write out a histogram of the number of bases in a capture covered at each depth or greater.
    > -Q:    Only count bases of at least <parameter> quality [0]
    > -f:    Only count bases of at least <parameter> Family size (unmarked reads are treated as FM 1) [0]
    > -m:    Max depth. Default: 262144.
    > -n:    Set N for quantile reporting. Default: 4 (quartiles)
    > -p:    Number of bases around region to pad in coverage calculations. Default: 0
    > -s:    Skip reads with an FP tag whose value is 0. (Fail)


####bmftools depth
Calculates depth of coverage across a bed file using barcode metadata.

####bmftools target
Calculates on-target fraction for bed file using barcode metadata.

####bmftools err
Calculates error rates by a variety of parameters.
Additionally, pre-computes the quality score recalibration for the optional dmp/sdmp recalibration step.

####bmftools famstats
Calculates summary statistics related to family size and demultiplexing.

###"Undocumented" tools

####<b>inmem</b>
  Description:
  > Largely an in-memory clone of bmftools dmp. Currently, only paired-end inline chemistry is supported.
  > In addition, use of rescaled quality scores is not yet supported.
  > RAM-hungry but fast. Useful for relatively small datasets.


  Usage: `bmftools inmem <options> input_R1.fast1.gz input_R2.fastq.gz`

  Options:

    > -1:    Output read1 path. REQUIRED.
    > -2:    Output read2 path. REQUIRED.
    > -=:    Emit output to stdout, interleaved if paired-end, instead of writing to disk.
    > -s:    Homing sequence. REQUIRED.
    > -l:    Barcode length. REQUIRED. For variable-length barcodes, this signifies the minimum length.
    > -v:    Maximum barcode length. Only needed for variable-length barcodes.
    > -t:    Reads with a homopolymer of threshold <parameter> length or greater are marked as QC fail. Default: 10.
    > -m:    Skip first <parameter> bases at the beginning of each read for use in barcode due to their high error rates.
    > -h/-?: Print help menu.

### Workflow
