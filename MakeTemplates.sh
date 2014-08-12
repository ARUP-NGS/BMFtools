#!/bin/bash

fq1=$1
fq2=$2
sample=$3
template="/yggdrasil/workspace/JavaWS/Pipeline/SarcomaDNAxLocTemplatev0.xml"

cat $template | sed "s:INPUTFILE1:$fq1:g" | sed "s:INPUTFILE2:$fq2:g" | sed "s:SAMPLE:$3:g" | sed 's:${::g' | sed 's:}::g' > ${3%.*}.Sarcoma.xml
