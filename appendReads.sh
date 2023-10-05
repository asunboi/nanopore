#!/bin/bash

query=$1
readnames=$2
outfile=$3

## example appendReads.sh 1_filter.fasta 2_readnames.txt 

faSomeRecords $query $readnames temp.fasta -exclude
cat temp.fasta >> $3