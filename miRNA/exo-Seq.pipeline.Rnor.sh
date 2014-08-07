#!/bin/bash
#5.8.2012
#This is Daniel Baker redoing a pipeline with shell and for Rattus norvegicus for miRNA

#Steps

#1. quality trim to QV 17
#Assume that adapters have already been removed. If not, use cutadapt
#2. remove 'too short, <17 nt' reads
#3. align to miRBase precursors with bowtie2
#3a. count and summarize mapping stats
#4. align unmapped reads from step 4 to frnadb, count, and summarize
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

ens_toplevel_ref="/mounts/genome/Rnor/Rattus_norvegicus.Rnor_5.0.75.dna.toplevel.fa"
mir_ref="/mounts/genome/Rnor/miRbase_v19_Rattus_norvegicus.fasta"
refseq_ref="/mounts/genome/Rnor/rat.rna.fna"
frna_ref="/mounts/genome/Rnor/Rattus_norvegicus.frnadb.fasta"
    #frna contains piRNA, miRNA, snoRNA, scRNA,snRNA,scaRNA,ribosomal, and other
    #No noncode entry for Rattus norvegicus
ens_ref="/mounts/genome/Rnor/Rattus_norvegicus.Rnor_5.0.75.ncrna.rename.fa"
#$noncode_ref="/yggdrasil/workspace/miRNA/miRNA_Analysis/miRNA_data/uniqncrna.2.0.human.fa"

#Begin alignments
$bt2 $mir_ref -U ../$fq --un $1.unmapped.mir.fastq -p $threads -S $1.mapped.mir.sam &> $1.mir.log

#Count hits
grep -v '@SQ\|@HD\|@PG' $1.mapped.mir.sam | grep 'AS:i' | awk '{print $3}' > $1.mapped.mir.names
cat $1.mapped.mir.names | sort | uniq -c | sort -k1,1n - | awk 'BEGIN {OFS="\t"};{print $1,$2}' > $1.mapped.mir.counts
echo 'total reads mapped to miRBase' $(wc -l $1.mapped.mir.names) > $1.report
echo 'number of genes/reference sequences detection for miRBase: ' $(wc -l $1.mapped.mir.counts) >> $1.report
echo 'top ten genes/references sequences in miRBase: ' >> $1.report
tail -n 10 $1.mapped.mir.counts >> $1.report

#Align to fRNAdb
$bt2 $frna_ref -U $fq.unmapped.mir.fastq --un $1.unmapped.mir.frna.fastq -p $threads -S $1.mapped.frna.sam &> $1.frna.log

#Count hits
grep -v '@SQ\|@HD\|@PG' $1.mapped.frna.sam | grep 'AS:i' | awk '{print $3}' > $1.mapped.frna.names
cat $1.mapped.frna.names | sort | uniq -c | sort -k1,1n -| awk 'BEGIN {OFS="\t"};{print $1,$2}' > $1.mapped.frna.counts
echo 'total reads mapped to fRNAdb' $(wc -l $1.mapped.frna.names) >> $1.report
echo 'number of genes/reference sequences detection for fRNAdb: ' $(wc -l $1.mapped.frna.counts) >> $1.report
echo 'top ten genes/references sequences in fRNAdb: ' >> $1.report
tail -n 10 $1.mapped.frna.counts >> $1.report

#Align to ensembl
$bt2 $ens_ref -U $fq.unmapped.mir.frna.fastq --un $1.unmapped.mir.frna.ens.fastq -p $threads -S $1.mapped.ens.sam &> $1.ens.log

#Count hits
grep -v '@SQ\|@HD\|@PG' $1.mapped.ens.sam | grep 'AS:i' | awk '{print $3}' > $1.mapped.ens.names
cat $1.mapped.ens.names | sort | uniq -c | sort -k1,1n - | awk 'BEGIN {OFS="\t"};{print $1,$2}' > $1.mapped.ens.counts
echo 'total reads mapped to ensembl' $(wc -l $1.mapped.ens.names) >> $1.report
echo 'number of genes/reference sequences detection for ensemble: ' $(wc -l $1.mapped.ens.counts) >> $1.report
echo 'top ten genes/references sequences in ensembl: ' >> $1.report
tail -n 10 $1.mapped.ens.counts >> $1.report

echo 'total reads unmapped to miRNA/fRNA/ensembl small RNA databases is '$(wc -l $fq.unmapped.mir.frna.ens.fastq) >> $1.report

#Align to RefSeq
$bt2 $refseq_ref -U $fq.unmapped.mir.frna.ens.fastq --un $1.unmapped.mir.frna.ens.rs.fastq -p $threads -S $1.mapped.rs.sam &> $1.rs.log

#Count hits
grep -v '@SQ\|@HD\|@PG' $1.mapped.rs.sam | grep 'AS:i' | awk '{print $3}' > $1.mapped.rs.names
cat $1.mapped.rs.names | sort | uniq -c | sort -k1,1n - | awk 'BEGIN {OFS="\t"};{print $1,$2}' > $1.mapped.rs.counts
echo 'total reads mapped to RefSeq' $(wc -l $1.mapped.rs.names) >> $1.report
echo 'number of genes/reference sequences detection for RefSeq: ' $(wc -l $1.mapped.rs.counts) >> $1.report
echo 'top ten genes/references sequences in RefSeq: ' >> $1.report
tail -n 10 $1.mapped.rs.counts >> $1.report

un=$(wc -l < $fq.unmapped.mir.frna.ens.rs.fastq)
und=$(expr $un / 4)
echo 'total reads completed unaligned: ' $und >> $1.report
tot=$(wc -l < ../$1)
totd=$(expr $tot / 4)
echo 'of total ' $(wc -l ../$1) >> $1.report
