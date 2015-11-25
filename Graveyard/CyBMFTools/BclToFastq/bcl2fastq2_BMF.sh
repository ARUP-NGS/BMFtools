#!/bin/bash

#Sample sheet must be prepared with the sample barcodes included.

OUTPUT_DIR=$1

/usr/local/bin/bcl2fastq2 -o $1 -R ./ --ignore-missing-bcls --ignore-missing-filter --ignore-missing-locs --use-bases-mask Y*,Y*,I8,Y* --minimum-trimmed-read-length 16

