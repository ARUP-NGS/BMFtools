#!/bin/bash
PT="/mounts/bin/picard-tools"
US="/yggdrasil/codeRepository/Packages/USeq_8.8.0/Apps"
Input=$1
InputPrefix=${Input%%.*}
InputBase=$(basename $InputPrefix)
bedfile="/mounts/Project_Sarcoma_Lib_709/sarcoma.bed"
SampDir=${Input%/*}
#SVDif=$2/$InputPrefix

java -jar $PT/MergeSamFiles.jar I=$InputPrefix.full/passSingle_$InputBase.full.ds.bam I=$InputPrefix.full/passSpan_$InputBase.full.ds.bam I=$InputPrefix.full/passSoft_$InputBase.full.ds.bam I=$InputPrefix.N30/passSpan_$InputBase.N30.ds.bam I=$InputPrefix.N30/passSingle_$InputBase.N30.ds.bam MSD=true AS=true O=$SampDir/merge$InputBase.bam
samtools index $SampDir/merge$InputBase.bam
delly -t 'TRA' -o $SampDir/merge$InputBase.vcf -g /mounts/genome/human_g1k_v37.fasta $SampDir/merge$InputBase.bam 
