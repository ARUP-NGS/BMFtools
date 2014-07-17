#!/bin/bash
BASENAME=${1%%.*}
cat ${BASENAME}.match.12.tags.fastq | paste - - - - | awk 'BEGIN {FS="\t";OFS="\t"};{print $2}' | sort | uniq -c | awk 'BEGIN {OFS="\t"};{print $1,$2}' > ${BASENAME}.fams
