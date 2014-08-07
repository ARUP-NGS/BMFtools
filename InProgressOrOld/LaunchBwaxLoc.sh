#!/bin/bash
#ref="/mnt/research1/reference/human_g1k_v37.fasta"
ref="/mounts/genome/human_g1k_v37.fasta"
threads=2
n=30
#PT=/mnt/research2/Daniel/bin/jar/picard-tools-1.109/picard-tools-1.109/
#US=/mnt/research2/Daniel/bin/jar/USeq_8.7.8/Apps/
bedfile="/mounts/Project_Sarcoma_Lib_709/sarcoma.bed"
Input1=$1
Input2=$2
InputPrefix=${Input1%.*}
InputAP=${Input1%.*}
bwa mem $ref $1 $2 -M -t $threads > $InputPrefix.full.sam
#Run full alignment

fastx_trimmer -i $1 -l $n -o $InputPrefix.N$n.fastq -Q33
fastx_trimmer -i $2 -l $n -o $InputAP.N$n.fastq -Q33
bwa mem $ref $InputPrefix.N$n.fastq $InputAP.N$n.fastq -t $threads > $InputPrefix.N$n.sam
$Run trimmed alignment

samtools view -Sbh $InputPrefix.full.sam > $InputPrefix.full.bam
samtools-mt sort -@ 2 $InputPrefix.full.bam  $InputPrefix.full.cso
java -Xmx3G -jar $PT/MarkDuplicates.jar I=$InputPrefix.full.cso.bam O=$InputPrefix.full.dd.bam REMOVE_DUPLICATES=true M=$InputPrefix.full.dd.log AS=true
samtools-mt sort -@ 2 -n $InputPrefix.full.dd.bam $InputPrefix.full.ds
java -Xmx8G -jar $US/SamSVFilter -a $InputPrefix.full.ds.bam -b $bedfile -s $InputPrefix.full

samtools view -Sbh $InputPrefix.N$n.sam > $InputPrefix.N$n.bam
samtools-mt sort -@ 2 $InputPrefix.N$n.bam  $InputPrefix.N$n.cso
java -Xmx3G -jar $PT/MarkDuplicates.jar I=$InputPrefix.N$n.cso.bam O=$InputPrefix.N$n.dd.bam REMOVE_DUPLICATES=true M=$InputPrefix.N$n.dd.log AS=true
samtools-mt sort -@ 2 -n $InputPrefix.N$n.dd.bam $InputPrefix.N$n.ds
java -Xmx8G -jar $US/SamSVFilter -a $InputPrefix.N$n.ds.bam -b $bedfile -s $InputPrefix.N$n

InputBase=$(basename $InputPrefix)
SampDir=${Input1%/*}

java -jar $PT/MergeSamFiles.jar I=$InputPrefix.N$n/passSingle$InputBase.N$n.ds.bam I=$InputPrefix.N$n/passSpan$InputBase.N$n.ds.bam I=$InputPrefix.full/passSingle$InputBase.full.ds.bam I=$InputPrefix.full/passSoft$InputBase.full.ds.bam I=$InputPrefix.full/passSpan$InputBase.full.ds.bam MSD=true AS=true USE_THREADING=true O=$SampDir/merge$InputBase.bam
samtools index $SampDir/merge$InputBase.bam

delly -t 'TRA' -g $ref -o $SampDir/merge$InputBase.vcf $SampDir/merge$InputBase.bam
