#!/bin/bash
PT="/mounts/bin/picard-tools"
US="/yggdrasil/codeRepository/Packages/USeq_8.8.0/Apps"
Input=$1
InputPrefix=${Input%.*}
bedfile="/mounts/Project_Sarcoma_Lib_709/sarcoma.bed"
samtools view -Sbh $1 > $InputPrefix.bam
samtools-mt sort -@ 2 $InputPrefix.bam  $InputPrefix.cso
java -Xmx3G -jar $PT/MarkDuplicates.jar I=$InputPrefix.cso.bam O=$InputPrefix.dd.bam REMOVE_DUPLICATES=true M=$InputPrefix.dd.log AS=true
samtools-mt sort -@ 2 -n $InputPrefix.dd.bam $InputPrefix.ds 
java -Xmx8G -jar $US/SamSVFilter -a $InputPrefix.ds.bam -b $bedfile -s $InputPrefix

rm $InputPrefix.cso.bam
rm $InputPrefix.dd.bam
