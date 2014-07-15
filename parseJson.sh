#!/bin/bash
FILENAME=$1
cat $FILENAME | tr ',' '\n' | sed 's:{:{\n:g' | sed 's:"\|{\|}::g' | sed 's/:/\t/g' | grep 'count' > ${FILENAME%%.*}.xTractCount
cat $FILENAME | tr ',' '\n' | sed 's:{:{\n:g' | sed 's:"\|{\|}::g' | sed 's/:/\t/g' | grep 'fraction' > ${FILENAME%%.*}.xTractFrac
