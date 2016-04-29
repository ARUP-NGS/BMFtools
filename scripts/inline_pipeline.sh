#!/bin/bash
set -e
set -x
REF="/mounts/genome/human_g1k_v37.fasta"
SORTMEM="4G"
SORT_THREADS1="2"
SORT_THREADS2="2"
THREADS="6"
BLEN=10
MAX_BLEN=13
r1=$1
r2=$2
HOMING="TGACT"
if [ -z "$PREFIX_LEN"]
then
    PREFIX_LEN="4"
fi
TMP_PREF="tmpfileswtf"
tmpstr=${1%.fq*}
FINAL_FQ_PREFIX=FINAL_${tmpstr%.fastq*}

# Perform inline barcode demultiplexing.
time bmftools dmp -zdp${THREADS} -l${BLEN} -v${MAX_BLEN} -s${HOMING} -n${PREFIX_LEN} \
	-o${TMP_PREF} $r1 $r2 -f${FINAL_FQ_PREFIX}
R1=${FINAL_FQ_PREFIX}.R1.fq.gz
R2=${FINAL_FQ_PREFIX}.R2.fq.gz
tmpstr=${R1%.fq*}
TMPFQ=${tmpstr%.fastq*}.tmp.fq
TMPBAM=${tmpstr%.fastq*}.tmp.bam
FINALBAM=${tmpstr%.fastq*}.rsq.bam

mkdir -p _FASTQC_${tmpstr}
mkdir -p _FASTQC_DMP_${tmpstr}
fastqc -t ${THREADS} --nogroup -o _FASTQC_${tmpstr} $1 $2
fastqc -t ${THREADS} --nogroup -o _FASTQC_DMP_${tmpstr} $R1 $R2

# There are a lot of processes here. We save a lot of time by avoiding I/O by piping.
bwa mem -CYT0 -t${THREADS} $REF $R1 $R2 | \
    bmftools mark -l 0 | \
    bmftools sort -l 0 -m $SORTMEM -@ $SORT_THREADS1 -k ucs -T tmpfileswtf | \
    bmftools rsq -usf $TMPFQ -l 0 - - | \
    samtools sort -O bam -T tmplastsort -@ $SORT_THREADS2 -m $SORTMEM -o $TMPBAM -

# Align the records that were rescued and merge them back in.
bwa mem -pCYT0 -t${THREADS} $REF $TMPFQ | bmftools mark -l 0 | \
    samtools sort -l 0 -O bam -T tmprsqsort -O bam -@ $SORT_THREADS2 -m $SORTMEM - | \
    samtools merge -fh $TMPBAM $FINALBAM $TMPBAM -

samtools index $FINALBAM

bmftools depth -sb $BED -p 50 $FINALBAM > ${FINALBAM}.doc.bed
bmftools famstats fm $TMPBAM > ${TMPBAM}.famstats.txt
bmftools famstats fm $FINALBAM > ${FINALBAM}.famstats.txt

# Clean up

rm $TMPFQ $TMPBAM
