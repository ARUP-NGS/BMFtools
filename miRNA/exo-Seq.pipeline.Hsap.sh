#!/bin/bash
#5.8.2012
#This is Daniel Baker redoing a pipeline with shell and for Rattus norvegicus for miRNA

#Steps

#1. quality trim to QV 17
#Assume that adapters have already been removed. If not, use cutadapt
#2. remove 'too short, <17 nt' reads
#3. align to miRBase precursors with bowtie2
#3a. count and summarize mapping stats
#4. align unmapped reads from step 4 to ENSncdb, count, and summarize
#5. align unmapped reads to ensembl, count, summarize
#6. align unmapped reads to RefSeq, count, summarize
#7. complete summary

fq=$1
outdir=$2
threads=4;

mkdir $outdir

echo "$outdir"
cd $outdir

bt2="/mounts/bin/bowtie2 --local --very-sensitive -p 14 --mm -x "
fastxtk="/mounts/bin/fastx_toolkit/src/fastq_quality_trimmer/fastq_quality_trimmer -v -t 17 -Q 33 -i "



mir_ref="/yggdrasil/workspace/miRNA/miRNA_Analysis/miRNA_data/miRBase.19.0.precursors.otherRNA.fasta";
refseq_ref="/yggdrasil/workspace/miRNA/miRNA_Analysis/miRNA_data/human.rna.fna";
noncode_ref="/yggdrasil/workspace/miRNA/miRNA_Analysis/miRNA_data/uniqncrna.2.0.human.fa";

#Noncoding Ensembl
ensnc_ref="/yggdrasil/workspace/miRNA/miRNA_Analysis/miRNA_data/Homo_sapiens.GRCh37.75.ncrna.fa"
#Coding Ensembl
ens_ref="/yggdrasil/workspace/miRNA/miRNA_Analysis/miRNA_data/Homo_sapiens.GRCh37.75.cds.all.fa"

enst_table="/yggdrasil/workspace/miRNA/miRNA_Analysis/miRNA_data/Ens2Gene.list"
enst_script="/yggdrasil/workspace/UnexpectedProtocols//ENST2Gene.py"



#Begin alignments
$bt2 $mir_ref -U ../$fq --un $1.unmapped.mir.fastq -p $threads -S $1.mapped.mir.sam &> $1.mir.log

#Count hits
grep -v '@SQ\|@HD\|@PG' $1.mapped.mir.sam | grep 'AS:i' | awk '{print $3}' > $1.mapped.mir.names
cat $1.mapped.mir.names | sort | uniq -c | sort -k1,1n - | awk 'BEGIN {OFS="\t"};{print $1,$2}' > $1.mapped.mir.counts
echo 'total reads mapped to miRBase' $(wc -l $1.mapped.mir.names) > $1.report
echo 'number of genes/reference sequences detection for miRBase: ' $(wc -l $1.mapped.mir.counts) >> $1.report
echo 'top ten genes/references sequences in miRBase: ' >> $1.report
tail -n 10 $1.mapped.mir.counts >> $1.report

#Align to Ensembl Non-coding ref
$bt2 $ensnc_ref -U $fq.unmapped.mir.fastq --un $1.unmapped.mir.ENSnc.fastq -p $threads -S $1.mapped.ENSnc.sam &> $1.ENSnc.log

#Count hits
grep -v '@SQ\|@HD\|@PG' $1.mapped.ENSnc.sam | grep 'AS:i' | awk '{print $3}' > $1.mapped.ENSnc.names
cat $1.mapped.ENSnc.names | sort | uniq -c | sort -k1,1n - | awk 'BEGIN {OFS="\t"};{print $1,$2}' > $1.mapped.ENSnc.counts
echo 'total reads mapped to ENSncdb' $(wc -l $1.mapped.ENSnc.names) >> $1.report
echo 'number of genes/reference sequences detection for ENSncdb: ' $(wc -l $1.mapped.ENSnc.counts) >> $1.report
echo 'top ten genes/references sequences in ENSncdb: ' >> $1.report
tail -n 10 $1.mapped.ENSnc.counts >> $1.report

#Align to uniqncrna ref
$bt2 $ensnc_ref -U $fq.unmapped.mir.ENSnc.fastq --un $1.unmapped.mir.ENSnc.uniqnc.fastq -p $threads -S $1.mapped.uniqnc.sam &> $1.uniqnc.log

#Count hits
grep -v '@SQ\|@HD\|@PG' $1.mapped.uniqnc.sam | grep 'AS:i' | awk '{print $3}' > $1.mapped.uniqnc.names
cat $1.mapped.uniqnc.names | sort | uniq -c | sort -k1,1n - | awk 'BEGIN {OFS="\t"};{print $1,$2}' > $1.mapped.uniqnc.counts
echo 'total reads mapped to NONCODE db' $(wc -l $1.mapped.uniqnc.names) >> $1.report
echo 'number of genes/reference sequences detection for NONCODE db: ' $(wc -l $1.mapped.uniqnc.counts) >> $1.report
echo 'top ten genes/references sequences in NONCODE db: ' >> $1.report
tail -n 10 $1.mapped.uniqnc.counts >> $1.report

#Align to ensembl
$bt2 $ens_ref -U $fq.unmapped.mir.ENSnc.uniqnc.fastq --un $1.unmapped.mir.ENSnc.uniqnc.ens.fastq -p $threads -S $1.mapped.ens.sam &> $1.ens.log

#Count hits
grep -v '@SQ\|@HD\|@PG' $1.mapped.ens.sam | grep 'AS:i' | awk '{print $3}' > $1.mapped.ens.names
cat $1.mapped.ens.names | sort | uniq -c | sort -k1,1n - | awk 'BEGIN {OFS="\t"};{print $1,$2}' > $1.mapped.ens.counts
echo 'total reads mapped to ensembl' $(wc -l $1.mapped.ens.names) >> $1.report
echo 'number of genes/reference sequences detection for ensemble: ' $(wc -l $1.mapped.ens.counts) >> $1.report
echo 'top ten genes/references sequences in ensembl: ' >> $1.report
tail -n 10 $1.mapped.ens.counts >> $1.report

echo 'total reads unmapped to miRNA/ENSnc/ensembl small RNA databases is '$(wc -l $fq.unmapped.mir.ENSnc.uniqnc.ens.fastq) >> $1.report

#Align to RefSeq
$bt2 $refseq_ref -U $fq.unmapped.mir.ENSnc.uniqnc.ens.fastq --un $1.unmapped.mir.ENSnc.uniqnc.ens.rs.fastq -p $threads -S $1.mapped.rs.sam &> $1.rs.log

#Count hits
grep -v '@SQ\|@HD\|@PG' $1.mapped.rs.sam | grep 'AS:i' | awk '{print $3}' > $1.mapped.rs.names
cat $1.mapped.rs.names | sort | uniq -c | sort -k1,1n - | awk 'BEGIN {OFS="\t"};{print $1,$2}' > $1.mapped.rs.counts
echo 'total reads mapped to RefSeq' $(wc -l $1.mapped.rs.names) >> $1.report
echo 'number of genes/reference sequences detection for RefSeq: ' $(wc -l $1.mapped.rs.counts) >> $1.report
echo 'top ten genes/references sequences in RefSeq: ' >> $1.report
tail -n 10 $1.mapped.rs.counts >> $1.report

un=$(wc -l < $fq.unmapped.mir.ENSnc.uniqnc.ens.rs.fastq)
und=$(expr $un / 4)
echo 'total reads completed unaligned: ' $und >> $1.report
tot=$(wc -l < ../$1)
totd=$(expr $tot / 4)
echo 'of total ' $(wc -l ../$1) >> $1.report

basename=$fq

cat ${basename}.mapped.ENSnc.counts | awk 'BEGIN {OFS="\t"};{print $0,"ENSncDB"}' > ${basename}.tmp.ENSnc
mv ${basename}.tmp.ENSnc ${basename}.mapped.ENSnc.counts
cat ${basename}.mapped.ens.counts | awk 'BEGIN {OFS="\t"};{print $0,"Ensembl"}' > ${basename}.tmp.ens
mv ${basename}.tmp.ens ${basename}.mapped.ens.counts
cat ${basename}.mapped.rs.counts | awk 'BEGIN {OFS="\t"};{print $0,"RefSeq"}' > ${basename}.tmp.rs
mv ${basename}.tmp.rs ${basename}.mapped.rs.counts
cat ${basename}.mapped.mir.counts | awk 'BEGIN {OFS="\t"};{print $0,"miRBase"}' > ${basename}.tmp.mir
mv ${basename}.tmp.mir ${basename}.mapped.mir.counts
cat ${basename}.mapped.uniqnc.counts | awk 'BEGIN {OFS="\t"};{print $0,"NONCODE"}' > ${basename}.tmp.uniqnc
mv ${basename}.tmp.uniqnc ${basename}.mapped.uniqnc.counts

cat *counts | sort -k1,1n > ${basename}.merged.counts

python $enst_script -c ${basename}.merged.counts -t $enst_table
