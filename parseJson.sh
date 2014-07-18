#!/bin/bash
FILENAME=$1
echo $FILENAME:Count > ${FILENAME%.*}.xTractCount
echo $FILENAME:Fraction  > ${FILENAME%.*}.xTractFrac 
cat $FILENAME | tr ',' '\n' | sed 's:{:{\n:g' | sed 's:"\|{\|}::g' | sed 's/:/\t/g' | grep 'count' | awk 'BEGIN {FS="\t";OFS="\t"};{print $NF}' >> ${FILENAME%.*}.xTractCount
cat $FILENAME | tr ',' '\n' | sed 's:{:{\n:g' | sed 's:"\|{\|}::g' | sed 's/:/\t/g' | grep 'fraction' | awk 'BEGIN {FS="\t";OFS="\t"};{print $NF}' >> ${FILENAME%.*}.xTractFrac
cat $FILENAME | tr ',' '\n' | sed 's:{:{\n:g' | sed 's:"\|{\|}::g' | sed 's/:/\t/g' | grep 'count' | awk 'BEGIN {FS="\t";OFS="\t"};{print $1}' >> ${FILENAME%.*}.xTractCountLabels
cat $FILENAME | tr ',' '\n' | sed 's:{:{\n:g' | sed 's:"\|{\|}::g' | sed 's/:/\t/g' | grep 'fraction' | awk 'BEGIN {FS="\t";OFS="\t"};{print $1}' >> ${FILENAME%.*}.xTractFracLabels

