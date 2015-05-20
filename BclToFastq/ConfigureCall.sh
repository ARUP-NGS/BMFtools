#!/bin/bash
#Command for the shades protocol fastq generation:

printf "Welcome to dual index fastq read generation!\nIf this is your first step in the process, you may want to go make sure your sample sheet is correctly prepared.\n"

SampleName=$1
SampleSheet=$2
BarcodeLen=$3
if [[ -z $SampleSheet ]]
    then
    printf "SampleSheet not set - using default value of ./SampleSheet.csv. \nThis might be wrong, especially if you are doing a dual index run.\n"
    BarcodeLen=8;
fi
exit 0
if [[ -z $BarcodeLen ]]
    then
    echo "BarcodeLength not set - using default value of 8\n"
    BarcodeLen=8;
fi
exit 0

# For index size 8 and an i7 index. (Call from the directory containing the Data folder)
/illumina/software/CASAVA-1.8.0/bin/configureBclToFastq.pl --input-dir Data/Intensities/BaseCalls/ --output-dir $SampleName --sample-sheet NewSampleSheetMod.csv --fastq-cluster-count 1000000000 --ignore-missing-stats --ignore-missing-bcl --use-bases-mask Y148,N16,I8,Y* --force
# For generating a fastq which has the indices
/illumina/software/CASAVA-1.8.0/bin/configureBclToFastq.pl --input-dir Data/Intensities/BaseCalls/ --output-dir ${SampleName}_index --sample-sheet NewSampleSheetMod.csv --fastq-cluster-count 1000000000 --ignore-missing-stats --ignore-missing-bcl --use-bases-mask Y148,Y16,I8,N* --force

#For an i5 index, switch the middle two numbers. BAM
# For index size 8 and an i7 index. (Call from the directory containing the Data folder)
/illumina/software/CASAVA-1.8.0/bin/configureBclToFastq.pl --input-dir Data/Intensities/BaseCalls/ --output-dir $SampleName --sample-sheet NewSampleSheetMod.csv --fastq-cluster-count 1000000000 --ignore-missing-stats --ignore-missing-bcl --use-bases-mask Y148,Y16,I8,Y* --force
# For generating a fastq which has the indices
/illumina/software/CASAVA-1.8.0/bin/configureBclToFastq.pl --input-dir Data/Intensities/BaseCalls/ --output-dir ${SampleName}_index --sample-sheet NewSampleSheetMod.csv --fastq-cluster-count 1000000000 --ignore-missing-stats --ignore-missing-bcl --use-bases-mask Y148,Y16,I8,N* --force
