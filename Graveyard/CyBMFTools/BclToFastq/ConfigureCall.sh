#!/bin/bash
#Command for the shades protocol fastq generation:

printf "Welcome to dual index fastq read generation!\nIf this is your first step in the process, you may want to go make sure your sample sheet is correctly prepared.\n"

SampleName=$1
SampleSheet=$2
BarcodeLen=$3
NmerIndexPos=$4
if [[ -z $SampleName ]]
    then
    EW_UUID=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 16 | head -n 1)
    SampleName="DefaultSampleName.${EW_UUID}"
    printf "SampleName not set - using a default with a UUID: ${SampleName}.\n"
fi
if [[ -z $SampleSheet ]]
    then
    printf "SampleSheet not set - using default value of ./SampleSheet.csv. \nThis might be wrong, especially if you are doing a dual index run.\n"
    SampleSheet="./SampleSheet.csv";
fi
if [[ -z $BarcodeLen ]]
    then
    echo "BarcodeLength not set - using default value of 8\n"
    BarcodeLen=8;
fi
if [[ -z $NmerIndexPos ]]
    then
    echo "NmerIndexPos not set - assuming i7. Variable is i5 if set and treated as i7 if not.\n"
    BarcodeLen=8;

    #For an i5 index, switch the middle two numbers. BAM
    # For index size 8 and an i7 index. (Call from the directory containing the ./Data folder)
    /illumina/software/CASAVA-1.8.0/bin/configureBclToFastq.pl --input-dir ./Data/Intensities/BaseCalls/ --output-dir $SampleName --sample-sheet ${SampleSheet} --fastq-cluster-count 1000000000 --ignore-missing-stats --ignore-missing-bcl --use-bases-mask Y148,Y16,I${BarcodeLen},Y* --force
    # For generating a fastq which has the indices
    /illumina/software/CASAVA-1.8.0/bin/configureBclToFastq.pl --input-dir ./Data/Intensities/BaseCalls/ --output-dir ${SampleName}_index --sample-sheet ${SampleSheet} --fastq-cluster-count 1000000000 --ignore-missing-stats --ignore-missing-bcl --use-bases-mask Y148,Y16,I${BarcodeLen},N* --force
    echo "Done running - go get your fastqs!"
    exit 0
fi

echo "NmerIndexPos set - assuming i5. Variable is i5 if set and treated as i7 if not.\n"

# For index size 8 and an i7 index. (Call from the directory containing the Data folder)
/illumina/software/CASAVA-1.8.0/bin/configureBclToFastq.pl --input-dir ./Data/Intensities/BaseCalls/ --output-dir $SampleName --sample-sheet ${SampleSheet} --fastq-cluster-count 1000000000 --ignore-missing-stats --ignore-missing-bcl --use-bases-mask Y148,N16,I${BarcodeLen},Y* --force
# For generating a fastq which has the indices
/illumina/software/CASAVA-1.8.0/bin/configureBclToFastq.pl --input-dir ./Data/Intensities/BaseCalls/ --output-dir ${SampleName}_index --sample-sheet ${SampleSheet} --fastq-cluster-count 1000000000 --ignore-missing-stats --ignore-missing-bcl --use-bases-mask Y148,Y16,I${BarcodeLen},N* --force
echo "Done running - go get your fastqs!"
exit 0
