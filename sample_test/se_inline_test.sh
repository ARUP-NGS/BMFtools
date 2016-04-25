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
if [ =z "$1"]
then
    r1="input.fastq.gz"
else
    r1=$1
fi
HOMING="TGACT"
if [ -z "$PREFIX_LEN"]
then
    PREFIX_LEN="4"
fi
TMP_PREF="tmpfileswtf"
tmpstr=${r1%.fq*}
echo $tmpstr
FINAL_FQ_PREFIX=FINAL_${tmpstr%.fastq*}
R1=${FINAL_FQ_PREFIX}.R1.fq.gz
tmpstr=${R1%.fq*}
LOG=${tmpstr}.log
> ${LOG} # Empty the file
TMPFQ=${tmpstr%.fastq*}.tmp.fq
TMPBAM=${tmpstr%.fastq*}.tmp.bam
PRERSQBAM=${tmpstr%.fastq*}.prersq.bam
FINALBAM=${tmpstr%.fastq*}.rsq.bam

# Perform inline barcode demultiplexing.
time bmftools_db dmp -Szdp${THREADS} -l${BLEN} -v${MAX_BLEN} -s${HOMING} -n${PREFIX_LEN} \
	-o${TMP_PREF} -f${FINAL_FQ_PREFIX} $r1
echo Number of reads before dmp: $(zgrep -c '^+$' $r1) >> $LOG

echo FM sum after dmp: $(zcat $R1 | paste -d'~' - - - - | cut -f1 -d'~'  | cut -d':' -f9 | cut -f1 | paste -sd+ | bc -l) >> $LOG
echo RV sum after dmp: $(zcat $R1 | paste -d'~' - - - - | cut -f1 -d'~' | cut -f5 | cut -d':' -f3 | paste -sd+ | bc) >> $LOG
echo Number of collapsed observations after dmp: $(zgrep -c '^+$' $R1) >> $LOG

# There are a lot of processes here. We save a lot of time by avoiding I/O by piping.
bwa mem -CYT0 -t${THREADS} $REF $R1 | samtools view -Sbh - > ${R1}.postalign.bam
getsums.py ${R1}.postalign.bam 2>>$LOG
echo Post-mark read count -cF2816: $(samtools view -cF2816 ${R1}.postalign.bam) >> $LOG &
echo Post-mark read count -cF2304: $(samtools view -cF2304 ${R1}.postalign.bam) >> $LOG &
bmftools_db mark -Sl6 ${R1}.postalign.bam ${R1}.postmark.bam
../fp512_test ${R1}.postmark.bam omg.bam && rm -f omg.bam
echo bmftools_db sort -l 9 -m 3G -@ 10 -Sk ucs -T tmpfileswtf ${R1}.postmark.bam > $PRERSQBAM
bmftools_db sort -l 9 -m 3G -@ 10 -Sk ucs -T tmpfileswtf ${R1}.postmark.bam > $PRERSQBAM
getsums.py $PRERSQBAM 2>>$LOG &
echo Post-BMF-sort read count -cF2816: $(samtools view -cF2816 $PRERSQBAM) >> $LOG &
echo Post-BMF-sort read count -cF2304: $(samtools view -cF2304 $PRERSQBAM) >> $LOG &

echo bmftools_db rsq -Suf $TMPFQ -l 0 $PRERSQBAM - (Then samtools sort)
bmftools_db rsq -Suf $TMPFQ -l 0 $PRERSQBAM - | \
samtools sort -Obam -T tmplastsort -@ $SORT_THREADS2 -m $SORTMEM -o $TMPBAM -

echo Post-rescue, \# of fastq records to be realigned: $(grep -c '^+$' $TMPFQ) >> $LOG &
echo Post-rescue, before merge read count -cF2816: $(samtools view -cF2816 $TMPBAM) >> $LOG &
echo Post-rescue, before merge read count -cF2304: $(samtools view -cF2304 $TMPBAM) >> $LOG &

# Align the records that were rescued and merge them back in.
bwa mem -CYT0 -t${THREADS} $REF $TMPFQ | bmftools_db mark -l 0 | \
    samtools sort -l 0 -O bam -T tmprsqsort -O bam -@ $SORT_THREADS2 -m $SORTMEM - | \
    samtools merge -fh $TMPBAM $FINALBAM $TMPBAM -
getsums.py $TMPBAM 2>>$LOG
echo Post-rescue, merged-reads count -cF2816: $(samtools view -cF2816 $FINALBAM) >> $LOG &
echo Post-rescue, merged-reads count -cF2304: $(samtools view -cF2304 $FINALBAM) >> $LOG &

getsums.py $FINALBAM 2>>$LOG

samtools index $FINALBAM

# QC
bmftools_db depth -sb $BED -p 50 $FINALBAM > ${FINALBAM}.doc.bed
bmftools_db famstats fm $PRERSQBAM > ${PRERSQBAM}.famstats.txt
bmftools_db famstats fm $FINALBAM > ${FINALBAM}.famstats.txt

# Clean up

rm $TMPFQ $TMPBAM
