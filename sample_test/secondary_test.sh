#!/bin/bash
# Usage:
# ./argv[0] r1.fq r2.fq r3.fq # r3 is the real r2, r2 is index
set -e
set -x
REF="/mounts/genome/human_g1k_v37.fasta"
BED="/mnt/d3/ctDNA_V3_hotspot_SNVs_ROI_with_Colo829_dopes.bed"
SORTMEM="4G"
SORT_THREADS1="2"
SORT_THREADS2="2"
THREADS="10"
SALT=6
r1=$1
rindex=$2
r2=$3
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

# Perform inline barcode demultiplexing.
echo time bmftools_db sdmp -zdp${THREADS} -s${SALT} -n${PREFIX_LEN} -i $rindex -o${TMP_PREF} $r1 $r2 -f${FINAL_FQ_PREFIX}
time bmftools_db sdmp -zdp${THREADS} -s${SALT} -n${PREFIX_LEN} -i $rindex -o${TMP_PREF} $r1 $r2 -f${FINAL_FQ_PREFIX}
echo Number of reads before dmp: $(zgrep -c '^+$' $r1) >> $LOG

echo FM sum after dmp: $(zcat $R1 | paste -d'~' - - - - | cut -f1 -d'~'  | cut -d':' -f9 | cut -f1 | paste -sd+ | bc ) >> $LOG
echo RV sum after dmp: $(zcat $R1 | paste -d'~' - - - - | cut -f1 -d'~' | cut -f5 | cut -d':' -f3 | paste -sd+ | bc) >> $LOG
echo Number of collapsed observations after dmp: $(zgrep -c '^+$' $R1) >> $LOG

# There are a lot of processes here. We save a lot of time by avoiding I/O by piping.
# echo bwa mem -CYT0 -t${THREADS} $REF $R1 $R2 | samtools view -Sbh -o ${R1}.premark.bam -
echo "bwa mem -CYT0 -t${THREADS} $REF $R1 $R2 | samtools view -Sbho ${R1}.premark.bam -"
bwa mem -CYT0 -t${THREADS} $REF $R1 $R2 | samtools view -Sbho ${R1}.premark.bam -
#bwa mem -CYT0 -t${THREADS} $REF $R1 $R2 | \
#    samtools view -Sbho ${R1}.premark.bam -
echo getsums.py ${R1}.premark.bam >>$LOG
getsums.py ${R1}.premark.bam >>$LOG
echo Pre-mark read count -cF2816: $(samtools view -cF2816 ${R1}.premark.bam) >> $LOG &
echo Pre-mark read count -cF2304: $(samtools view -cF2304 ${R1}.premark.bam) >> $LOG &
echo bmftools_db mark ${R1}.premark.bam ${R1}.postmark.bam
bmftools_db mark ${R1}.premark.bam ${R1}.postmark.bam
getsums.py ${R1}.postmark.bam >>$LOG
echo Post-mark read count -cF2816: $(samtools view -cF2816 ${R1}.postmark.bam) >> $LOG &
echo Post-mark read count -cF2304: $(samtools view -cF2304 ${R1}.postmark.bam) >> $LOG &
bwa mem -CYT0 -t${THREADS} $REF $R1 $R2 | samtools view -Sbh - > ${R1}.premark.bam
echo getsums.py ${R1}.premark.bam >>$LOG &
getsums.py ${R1}.premark.bam >>$LOG &
echo Pre-mark read count -cF2816: $(samtools view -cF2816 ${R1}.premark.bam) >> $LOG &
echo Pre-mark read count -cF2304: $(samtools view -cF2304 ${R1}.premark.bam) >> $LOG &
echo bmftools_db mark ${R1}.premark.bam ${R1}.postmark.bam
bmftools_db mark ${R1}.premark.bam ${R1}.postmark.bam
getsums.py ${R1}.postmark.bam >>$LOG
echo Post-mark read count -cF2816: $(samtools view -cF2816 ${R1}.postmark.bam) >> $LOG &
echo Post-mark read count -cF2304: $(samtools view -cF2304 ${R1}.postmark.bam) >> $LOG &
echo "Now checing that the qcfail bit is set correctly."
echo fp512_test ${R1}.postmark.bam /dev/null
fp512_test ${R1}.postmark.bam /dev/null
echo bmftools_db sort -l 9 -m 3G -@ 10 -k ucs -T tmpfileswtf ${R1}.postmark.bam > $PRERSQBAM
bmftools_db sort -l 9 -m 3G -@ 10 -k ucs -T tmpfileswtf ${R1}.postmark.bam > $PRERSQBAM
getsums.py $PRERSQBAM >>$LOG &
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
getsums.py $TMPBAM >>$LOG
echo Post-rescue, merged-reads count -cF2816: $(samtools view -cF2816 $FINALBAM) >> $LOG &
echo Post-rescue, merged-reads count -cF2304: $(samtools view -cF2304 $FINALBAM) >> $LOG &

getsums.py $FINALBAM >>$LOG

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
