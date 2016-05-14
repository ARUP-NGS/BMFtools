#BMFtools
###Summary
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

###Synopsis
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

###Commands and options

#####<b>dmp</b>
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


#####<b>dmp</b>
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

### Workflow
