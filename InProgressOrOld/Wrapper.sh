#!/bin/bash
Full=$1
N30=$2

/mounts/Project_Sarcoma_Lib_709/BamUtils.sh $1
/mounts/Project_Sarcoma_Lib_709/BamUtils.sh $2
/mounts/Project_Sarcoma_Lib_709/MergeDelly.sh $1
