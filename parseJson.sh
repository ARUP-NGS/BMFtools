#!/bin/bash
FILENAME=$1
echo $FILENAME:Count > ${FILENAME%.*}.xTractCount
echo $FILENAME:Fraction  > ${FILENAME%.*}.xTractFrac 
cat $FILENAME | tr ',' '\n' | sed 's:{:{\n:g' | sed 's:"\|{\|}::g' | sed 's/:/\t/g' | grep 'count' >> ${FILENAME%.*}.xTractCount
cat $FILENAME | tr ',' '\n' | sed 's:{:{\n:g' | sed 's:"\|{\|}::g' | sed 's/:/\t/g' | grep 'fraction' >> ${FILENAME%.*}.xTractFrac
