#BMFtools
##Summary
BMFtools (**B**arcoded **M**olecular **F**amilies tools) is a suite of tools for error correction and precise quantitation using reads with molecular barcodes.
Reads which are PCR duplicates of original template molecules are **molecularly** demultiplexed into single unique observations for each sequenced founding template molecule.
Accessory tools provide postprocessing, filtering, quality control, and summary statistics.

##Index

[Workflow](#workflow)

[Tools](#tools)

[Usage](#usage)

See also [Vignettes](https://github.com/ARUP-NGS/BMFtools/blob/dev/Vignettes.md) in a separate document for use cases.

###Workflow

Reads are originally collapsed into single observations per barcode using exact matching. Following this step, to correct for errors in barcode reading, positional information can be used to rescue reads into unique single observations.

The only difference between inline and secondary index chemistry workflows is the initial `bmftools dmp` call.

####Exact-Matching Fastq Consolidation
A typical paired-end exact fastq-stage molecular demultiplexing call:

> Inline

```bmftools dmp -s <homing_sequence> -l <barcode_length> -o <temporary_file_prefix> -p <threads> -f <final_output_prefix> <r1.fq.gz> <r2.fq.gz>```
> Secondary Index

```bmftools sdmp -o <temporary_file_prefix> -p <threads> -f <final_output_prefix> -i <index.fq.gz> <r1.fq.gz> <r2.fq.gz>```

Each of these produces `final_output_prefix.R1.fq` and `final_output_prefix.R2.fq`. (A .gz suffix is appended if the output is gzip compressed.)

Barcode metadata is written into fastq comments in SAM auxiliary tag format. These reads are then aligned with
bwa mem with the -C option, which appends the fastq comment to the end of the sam record. This trivially adds
tags to all alignments for each read.

####Alignment

`bwa mem -CYT0 -t<threads> <idx.base> final_output_prefix.R1.fq final_output_prefix.R2.fq | samtools view -bho final_output.bam`


####Rescue
Because errors occur in reading barcodes, this initial exact-matching step is not completely successful in grouping
reads from the same original template molecule. To account for this, an optional rescue protocol has been implemented.

This rescue takes place in two steps -- first, a sort which groups together based on alignment signature, and second,
collapsing reads sharing these signatures with similar barcodes into single observations.

"Alignment signatures" consist of a read and its mate's alignment information, if paired. These can be grouped by start position
or by unclipped start position. Unclipped start position is less sensitive to errors in the reads, whereas
a bam sorted by signature using position can still be indexed for traditional use. Unclipped start position comes at the computational cost of an additional sort but with potentially increased success in rescue.

Because reads need both their and their mates' alignment information, including read length, the preprocessing
`bmftools mark` is required prior to bmftools sort.


Because reads that have been modified may align elsewhere, these reads are all realigned. In addition, this regenerates
all of the secondary and supplementary alignments.


In the collapsing `bmftools rsq` step, supplementary and secondary reads are stripped to preserve balanced pairs.
If secondary and supplementary alignments are needed for other reads,these should be written to the temporary fastq for realignment using the -s option.

`bmftools mark -l0 final_output.bam | sort -k <ucs/bmf> -o <final_output_prefix.bmfsort.bam> -`

For position:

`bmftools rsq -f<tmp.fq> <final_output_prefix.bmfsort.bam> <final_output_prefix.tmprsq.bam>`

For unclipped start:

`bmftools rsq [-u <unclipped start only>] -f<tmp.fq> <final_output_prefix.bmfsort.bam> - | samtools sort -O bam -T<tmp_prefix> -ofinal_output_prefix.tmprsq.bam`


Realigned reads are then sorted and merged in with the other reads in the dataset.

`bwa mem -pCYT0 -t<threasd> <reference> -f<tmp.fq> | bmftools mark |  samtools sort -l 0 -Obam -T <tmp_prefix> | samtools merge -cpfh final_output_prefix.tmprsq.bam final_output_prefix.rsqmerged.bam final_output_prefix.tmprsq.bam -`

For efficiency, this can be heavily piped to reduce I/O and unnecessary compression/decompression.
At this point, final_output_prefix.tmprsq.bam contains supplementary and secondary alignments for all template molecules, and reads have been rescued using alignment information.

####Post-rescue: Now what?

If you're not interested in taking advantage of the barcode metadata (new p values, family sizes, &c.), this bam is ready for downstream analysis.

For more specific use cases, postprocessing steps (such as cap and filter) can be used to prepare a bam for use by BMF-agnostic tools, and variant calling can be performed or vetted with `bmftools stack` or `bmftools vet`.

For translocation detection, soft-clippings are often used as markers for potential events. To improve performance of such tools, (e.g., WHAM), if adapters sequences have been masked, successive masked bases at the ends of reads can be removed by [maskripper](https://github.com/noseatbelts/maskripper). For rescue to work properly, however, this should be performed only after rescue has been completed.

##Tools


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

`bmftools cap -f 0.8 -c 200 input.bam output.bam`


`bmftools depth -H coverage_uniformity.hist.txt -b capture.bed -Q1 input.bam > coverage.bed.txt`


`bmftools dmp -p 12 -f output_prefix -l <barcode_length> input_R1.fastq.gz input_r2.fastq.gz`


`bmftools err fm reference.fasta input.srt.bam > err_by_fm.txt`


`bmftools err main -c cycle_err.txt -n base_call_err.txt -g global_err.txt -o rescaled_qualities.txt reference.fasta input.srt.bam`


`bmftools err region -b regions.bed -o err_by_region.txt reference.fasta input.srt.bam`


`bmftools famstats fm input.bam > famstats.txt`


`bmftools famstats frac 2 input.bam > frac_fm_ge2.txt`


`bmftools filter -s2 input.bam output.bam`


`bmftools mark input.bam output.bam`


`bmftools rsq -u -f tmp.fastq input.bam tmp.out.bam`


`bmftools sdmp -p 12 -s 1 -o tmp_prefix -f output_prefix -i index.fastq.gz read1.fastq.gz read2.fastq.gz`


`bmftools sort -T tmp_prefix -k ucs input.bam output.bam`


`bmftools stack --ref reference.fasta -o output.vcf -b capture.bed --min-family-size 3 tumor.bam normal.bam`


`bmftools target -b capture.bed input.bam`


`bmftools vet -o output.vcf -b capture.bed --min-family-size 3 input.bcf input.bam`

##Usage

### Core Functionality

####<b>dmp</b>
  Description:
  > Performs molecular demutiplexing of inline barcoded fastq data.
  > The homing sequence is a sequence of bases marking the end of the random nucleotides
  > which make up the barcode. This is required to ensure that chemistry has worked as expected
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
    > -h/-?: Print usage.


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
    > -h/-?: Print usage.


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
    > -h/-?: Print usage.

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
    * BMF_PASS: 1 if variant passes, 0 otherwise.
    * ADP: Number of unique observations for allele.
    * ADPO: Number of overlapping read pair observations for allele.
    * ADPD: Number of duplex observations for allele.
    * ADPR: Number of original reversed observations for allele.
    * RVF: Fraction of RV observations supporting allele.
    * QSS: Q Score Sum supporting allele.
    * AMBIG: Number of ambiguous base calls at position.
    * SOMATIC_CALL: Boolean value for a somatic call for allele.
    * FR_FAILED: Number of observations failed by fraction of family members agreed on a base call per sample.
    * FM_FAILED: Number of observations failed for insufficient family size per sample.
    * FP_FAILED: Number of observations failed for failing barcode QC per sample.
    * AF_FAILED: Number of observations failed for insufficient aligned fraction.
    * MQ_FAILED: Number of observations failed for insufficient mapping quality.
    * IMPROPER: Number of observations failed for being in an improper pair.
    * OVERLAP: Number of overlapping read pairs at position.
    * AFR: Allele Fractions per allele, including reference.

####<b>vet</b>
    Description:
    > Curates variant calls from a bcf file and an associated, indexed bam file.
    > Overlapping reads in a pileup are treated as two observations and are merged, if agreed, using Fisher's method.
    > Discordant unambiguous base calls are masked. Ambiguous base calls are overridden by an unambiguous call from its mate.

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

    VCF Header Fields:
    * BMF_VET:  1 if a variant passes, 0 otherwise. 
    * BMF_UNIOBS: Number of unique observations supporting a variant at that position.
    * BMF_DUPLEX: Number of duplex observations supporting a variant at that position.
    * BMF_FAIL: NUmber of reads at position failing filters.
    * DUPLEX_DEPTH: Number of duplex reads at position passing filters.
    * DISC_OVERLAP: Number of read pairs at position with discordant base calls.
    * OVERLAP: Number of overlapping read pairs combined into single observations at position.

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


####<b>target</b>
  Description:
  > Calculates the fraction of on-target reads, both raw and consolidated.

  Usage: bmftools target <opts> <in.bam>

  Options:

    > -b:    Path to bed. REQUIRED.
    > -m:    Set minimum mapping quality for inclusion.
    > -p:    Set padding - number of bases around target region to consider as on-target. Default: 0.
    > -n:    Set notification interval - number of reads between logging statements. Default: 1000000.

####<b>err</b>
  Description:
  > Calculates error rates by a variety of parameters.
  > Additionally, pre-computes the quality score recalibration for the optional dmp/sdmp recalibration step.
  > err has 3 subcommands:
  1. main
    1. Primary output is recalibrated quality scores given a sequenced, aligned, sorted standard dataset. (e.g., PhiX)
    2. err main also produces error rates by cycle, base call, quality score, facilitating error analysis.
  2. fm
    1. err fm calculates error rates by family size.
  3. region
    1. err region calculates error rates by bed region.

  Usage: bmftools err main <opts> <reference.fasta> <in.csrt.bam>

  Options:

    > -o:    Path to output file. Set to '-' or 'stdout' to emit to stdout.
    > -a:    Set minimum mapping quality for inclusion.
    > -S:    Set minimum calculated PV tag value for inclusion.
    > -r:    Name of contig. If set, only reads aligned to this contig are considered
    > -3:    Path to write the 3d offset array in tabular format.
    > -f:    Path to write the full measured error rates in tabular format.
    > -n:    Path to write the cycle/nucleotide call error rates in tabular format.
    > -c:    Path to write the cycle error rates in tabular format.
    > -g:    Path to write the global error rates in tabular format.
    > -b:    Path to bed file for restricting analysis.
    > -m:    Minimum family size for inclusion. Default: 0.
    > -M:    Maximum family size for inclusion. Default: 2147483647.
    > -d:    Flag to only calculate error rates for duplex reads.
    > -D:    Flag to only calculate error rates for non-duplex reads.
    > -p:    Set padding for bed region. Default: 0.
    > -P:    Only include proper pairs.
    > -O:    Set minimum number of observations for imputing quality Default: 10000.
    > -h/-?  Print usage.

  Usage: bmftools err fm <opts> <reference.fasta> <in.csrt.bam>

  Options:

    > -o:    Path to output file. Set to '-' or 'stdout' to emit to stdout.
    > -h/-?: Print usage.
    > -S:    Set minimum calculated PV tag value for inclusion.
    > -a:    Set minimum mapping quality for inclusion.
    > -r:    Name of contig. If set, only reads aligned to this contig are considered
    > -b:    Path to bed file for restricting analysis.
    > -d:    Flag to only calculate error rates for duplex reads.
    > -p:    Set padding for bed region. Default: 0.
    > -P:    Only include proper pairs.
    > -F:    Require that the FP tag be present and nonzero.
    > -f:    Require that the fraction of family members agreed on a base be <parameter> or greater. Default: 0.0

  Usage: bmftools err region <opts> <reference.fasta> <in.csrt.bam>

  Options:

   > -b:    Path to bed file. REQUIRED.
   > -o:    Path to output file. Leave unset or set to '-' or 'stdout' to emit to stdout.
   > -a:    Set minimum mapping quality for inclusion.
   > -p:    Set padding for bed region. Default: 0.
   > -h/-?: Print usage.


####bmftools famstats
  Description:
  > Calculates summary statistics related to family size and demultiplexing.
  > famstats consists of two subcommands: fm and frac
  > famstats fm has 2 subcommands:
  1. fm
    2. famstats fm produces summary statistics and count distributions for family size, duplex/reverse reads, and read rescue statistics.
  2. frac
    1. famstats frac

  Usage: bmftools famstats fm <opts> <in.bam>

  Options:

    > -m:    Set minimum mapping quality. Default: 0.
    > -f:    Set minimum family size. Default: 0.

  Usage: bmftools famstats frac <opts> <minFM> <in.bam>

  Options:
    > -n:    Set notification interval. Default: 1000000.
    > -h/-?: Print usage.


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
    > -h/-?: Print usage.

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



### Utilities

####bmftools sort
  Description:
  > Sorts an alignment file in preparation for read consolidation using positional information.
  > Essentially a modification of samtools sort.

  Options:

    > -l INT       Set compression level, from 0 (uncompressed) to 9 (best)
    > -m INT       Set maximum memory per thread; suffix K/M/G recognized [768M]
    > -k           Sort key - pos for positional (samtools default), qname for query name, bmf for extended positional, ucs for using unclipped mate start/stop positions. Default: bmf comparison.
    > -o FILE      Write final output to FILE rather than standard output. If splitting, this is used as the prefix.
    > -O FORMAT    Write output as FORMAT ('sam'/'bam'/'cram') Default: bam.
    > -T PREFIX    Write temporary files to PREFIX.nnnn.bam. Default: 'MetasyntacticVariable')
    > -@ INT       Set number of sorting and compression threads [1]
    > -s           Flag to split the bam into a list of file handles.
    > -p           If splitting into a list of handles, this sets the file prefix.
    > -S           Flag to specify single-end. Needed for unclipped start compatibility.
    > -h/-?        Print usage.

####bmftools mark
  Description:
  > Marks a sets of template bam records with auxiliary tags for use in downstream tools.
  > Required for sort and rsq.
  > Intended primarily for piping. Default compression is therefore 0. Typical compression for writing to disk: 6.

  Usage: bmftools mark <opts> <input.namesrt.bam> <output.bam>

  Options:

    > -l:    Sets bam compression level. (Valid: 1-9). Default: 0.
    > -q:    Skip read pairs which fail.
    > -d:    Set bam compression level to default (6).
    > -i:    Skip read pairs whose insert size is less than <INT>.
    > -u:    Skip read pairs where both reads have a fraction of unambiguous base calls >= <parameter>
    > -S:    Use this for single-end marking. Only sets the QC fail bit for reads failing barcode QC.
    > Set input.namesrt.bam to '-' or 'stdin' to read from stdin.
    > Set output.bam to '-' or 'stdout' or omit to stdout.
    > Thus `bmftools mark` defaults to reading and writing from stdin and stdout, respectively, in paired-end mode.


###"Undocumented" tools

####<b>inmem</b>
  Description:
  > Largely an in-memory clone of bmftools dmp. Currently, only paired-end inline chemistry is supported.
  > In addition, use of rescaled quality scores is not yet supported.
  > RAM-hungry but fast. Useful for relatively small datasets.

  Usage: `bmftools inmem <options> input_R1.fastq.gz input_R2.fastq.gz`

  Options:

    > -1:    Output read1 path. REQUIRED.
    > -2:    Output read2 path. REQUIRED.
    > -=:    Emit output to stdout, interleaved if paired-end, instead of writing to disk.
    > -s:    Homing sequence. REQUIRED.
    > -l:    Barcode length. REQUIRED. For variable-length barcodes, this signifies the minimum length.
    > -v:    Maximum barcode length. Only needed for variable-length barcodes.
    > -t:    Reads with a homopolymer of threshold <parameter> length or greater are marked as QC fail. Default: 10.
    > -m:    Skip first <parameter> bases at the beginning of each read for use in barcode due to their high error rates.
    > -h/-?: Print usage.


####<b>hashdmp</b>
  Description:

   > Contains the hashmap-powered consolidation only. Its input is the preprocessed marked bams produced
   > by bmftools dmp and sdmp.
   > Molecularly demultiplexes marked temporary fastqs into final unique observation records.
   > bmftools hashdmp does so in one large hashmap. This may require huge amounts of memory.


  Usage: `bmftools hashdmp <opts> <input.marked.fq>`

  Options:

    > -s:    Perform secondary index consolidation rather than Loeb-like inline consolidation.
    > -o:    Write to <parameter> rather than stdout.
