#!/bin/bash

cat $1 | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n"
