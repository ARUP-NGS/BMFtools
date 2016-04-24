#!/bin/bash
set -e
set -x
REF="/mounts/genome/human_g1k_v37.fasta"
BED="/mnt/d3/ctDNA_V3_hotspot_SNVs_ROI_with_Colo829_dopes.bed"
SORTMEM="4G"
SORT_THREADS1="2"
SORT_THREADS2="2"
THREADS="10"
BLEN=13
MAX_BLEN=13
HOMING="TGACT"
if [ -z "$PREFIX_LEN"]
then
    PREFIX_LEN="4"
fi
TMP_PREF="tmpfileswtf"
tmpstr=${1%.fq*}
FINAL_FQ_PREFIX=FINAL_${tmpstr%.fastq*}
R1=${FINAL_FQ_PREFIX}.R1.fq.gz
R2=${FINAL_FQ_PREFIX}.R2.fq.gz
tmpstr=${R1%.fq*}
LOG=${tmpstr}.log
> ${LOG} # Empty the file
TMPFQ=${tmpstr%.fastq*}.tmp.fq
TMPBAM=${tmpstr%.fastq*}.tmp.bam
PRERSQBAM=${tmpstr%.fastq*}.prersq.bam
FINALBAM=${tmpstr%.fastq*}.rsq.bam

echo Pre-align number of dmp reads: $(zgrep -c '^+$' $1)
python -c "import pysam; print sum([int(a.comment.split()[3].split(\":\")[2]) for a in pysam.FastqFile(\"${1}\")])"
# There are a lot of processes here. We save a lot of time by avoiding I/O by piping.
bwa mem -CYT0 -t${THREADS} $REF $R1 $R2 | samtools view -Sbh - > ${R1}.premark.bam
getsums.py ${R1}.premark.bam 2>>$LOG &
echo Pre-mark read count -cF2816: $(samtools view -cF2816 ${R1}.premark.bam) >> $LOG &
echo Pre-mark read count -cF2304: $(samtools view -cF2304 ${R1}.premark.bam) >> $LOG &
echo bmftools_db mark ${R1}.premark.bam ${R1}.postmark.bam
bmftools_db mark ${R1}.premark.bam ${R1}.postmark.bam
getsums.py ${R1}.postmark.bam 2>>$LOG
echo Post-mark read count -cF2816: $(samtools view -cF2816 ${R1}.postmark.bam) >> $LOG &
echo Post-mark read count -cF2304: $(samtools view -cF2304 ${R1}.postmark.bam) >> $LOG &
echo "Now checing that the qcfail bit is set correctly."
echo ../fp512 ${R1}.postmark.bam omg.bam && rm -f omg.bam
../fp512 ${R1}.postmark.bam omg.bam && rm -f omg.bam
echo bmftools_db sort -l 9 -m 3G -@ 10 -k ucs -T tmpfileswtf ${R1}.postmark.bam > $PRERSQBAM
bmftools_db sort -l 9 -m 3G -@ 10 -k ucs -T tmpfileswtf ${R1}.postmark.bam > $PRERSQBAM
getsums.py $PRERSQBAM 2>>$LOG &
echo Post-BMF-sort read count -cF2816: $(samtools view -cF2816 $PRERSQBAM) >> $LOG &
echo Post-BMF-sort read count -cF2304: $(samtools view -cF2304 $PRERSQBAM) >> $LOG &

bmftools_db rsq -uf $TMPFQ -l 0 $PRERSQBAM - | \
samtools sort -O bam -T tmplastsort -@ $SORT_THREADS2 -m $SORTMEM -o $TMPBAM -

echo Post-rescue, \# of fastq records to be realigned: $(grep -c '^+$' $TMPFQ) >> $LOG &
echo Post-rescue, before merge read count -cF2816: $(samtools view -cF2816 $TMPBAM) >> $LOG &
echo Post-rescue, before merge read count -cF2304: $(samtools view -cF2304 $TMPBAM) >> $LOG &

# Align the records that were rescued and merge them back in.
bwa mem -pCYT0 -t${THREADS} $REF $TMPFQ | bmftools_db mark -l 0 | \
    samtools sort -l 0 -O bam -T tmprsqsort -O bam -@ $SORT_THREADS2 -m $SORTMEM - | \
    samtools merge -fh $TMPBAM $FINALBAM $TMPBAM -
getsums.py $TMPBAM 2>>$LOG
echo Post-rescue, merged-reads count -cF2816: $(samtools view -cF2816 $TMPBAM) >> $LOG &
echo Post-rescue, merged-reads count -cF2304: $(samtools view -cF2304 $TMPBAM) >> $LOG &

# Align the records that were rescued and merge them back in.
bwa mem -pCYT0 -t${THREADS} $REF $TMPFQ | bmftools_db mark -l 0 | \
    samtools sort -l 0 -O bam -T tmprsqsort -O bam -@ $SORT_THREADS2 -m $SORTMEM - | \
    samtools merge -fh $TMPBAM $FINALBAM $TMPBAM -

getsums.py $FINALBAM 2>>$LOG

samtools index $FINALBAM

# QC


samtools index $FINALBAM

# QC

#bmftools_db depth -sb $BED -p 50 $PRERSQBAM > ${PRERSQBAM}.doc.bed
bmftools_db depth -sb $BED -p 50 $FINALBAM > ${FINALBAM}.doc.bed
bmftools_db famstats fm $PRERSQBAM > ${PRERSQBAM}.famstats.txt
bmftools_db famstats fm $FINALBAM > ${FINALBAM}.famstats.txt

# Clean up

rm $TMPFQ $TMPBAM
