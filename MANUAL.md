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

### Workflow
