#!/bin/bash
set -e
set -x
R1=$1
R2=$2
tmpstr=${1%.fq*}
TMPFQ=${tmpstr%.fastq*}.tmp.fq
TMPBAM=${tmpstr%.fastq*}.tmp.bam
FINALBAM=${tmpstr%.fastq*}.rsq.bam
REF="/mounts/genome/human_g1k_v37.fasta"
SORTMEM="4G"
SORT_THREADS="4"
THREADS="10"

bwa mem -CYT0 -t${THREADS} $REF $R1 $R2 | \
    bmftools mark_unclipped -l 0 - - | \
    bmftools sort -l 0 -m $SORTMEM -@ $SORT_THREADS -k ucs -T tmpfileswtf | \
    bmftools rsq -uf $TMPFQ -l 0 - - | \
    samtools sort -O bam -T tmplastsort -@ $SORT_THREADS -m $SORTMEM -o $TMPBAM -

# Sort fastq by name, align, sort/convert, merge with temporary bam
cat $TMPFQ | paste -d'~' - - - - | sort | tr '~' '\n' | \
    bwa mem -pCYT0 -t${THREADS} $REF - | \
    samtools sort -l 0 -O bam -T tmprsqsort -O bam -@ $SORT_THREADS -m $SORTMEM - | \
    samtools merge -h $TMPBAM $FINALBAM $TMPBAM -

# Clean up

rm $TMPFQ $TMPBAM
