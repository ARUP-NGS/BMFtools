#!/bin/bash
set -e
set -x
REF="/mounts/genome/human_g1k_v37.fasta"
SORTMEM="4G"
SORT_THREADS1="1"
SORT_THREADS2="4"
THREADS="10"
BLEN=10
MAX_BLEN=11
r1=$1
r2=$2
HOMING="TGACT"
PREFIX_LEN="3"
TMP_PREF="tmpfileswtf"
tmpstr=${1%.fq*}
FINAL_FQ_PREFIX=FINAL_${tmpstr%.fastq*}

time bmftools dmp -zdp${THREADS} -l${BLEN} -v${MAX_BLEN} -s${HOMING} -n${PREFIX_LEN} \
	-o${TMP_PREF} $r1 $r2 -f${FINAL_FQ_PREFIX}
R1=${FINAL_FQ_PREFIX}.R1.fq.gz
R2=${FINAL_FQ_PREFIX}.R2.fq.gz
tmpstr=${R1%.fq*}
TMPFQ=${tmpstr%.fastq*}.tmp.fq
TMPBAM=${tmpstr%.fastq*}.tmp.bam
FINALBAM=${tmpstr%.fastq*}.rsq.bam

bwa mem -CYT0 -t${THREADS} $REF $R1 $R2 | \
    bmftools mark -l 0 - - | \
    bmftools sort -l 0 -m $SORTMEM -@ $SORT_THREADS1 -k ucs -T tmpfileswtf | \
    bmftools rsq -uf $TMPFQ -l 0 - - | \
    samtools sort -O bam -T tmplastsort -@ $SORT_THREADS2 -m $SORTMEM -o $TMPBAM -

# Sort fastq by name, align, sort/convert, merge with temporary bam
bwa mem -pCYT0 -t${THREADS} $REF $TMPFQ | \
    samtools sort -l 0 -O bam -T tmprsqsort -O bam -@ $SORT_THREADS2 -m $SORTMEM - | \
    samtools merge -h $TMPBAM $FINALBAM $TMPBAM -

# Clean up

rm $TMPFQ $TMPBAM
