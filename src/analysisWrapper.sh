#!/bin/bash

mkdir $1
mkdir $1/cleaning
mkdir $1/clustering

### check whether the compiled fasta has already been made - if not, make it.
dir_name=`pwd`
barcode=${dir_name##*/}

if ! test -f ${barcode}.fa; then
  cat *.gz > ${barcode}.fastq.gz
  gunzip ${barcode}.fastq.gz
  seqtk seq -a ${barcode}.fastq > ${barcode}.fa
fi

cd $1
sbatch ~/scripts/nanopore/analysis.sbatch $1
cd ..

if [ ! -d "tong" ]; then
    mkdir tong
    cd "tong" || exit 1
    sbatch ~/scripts/nanopore/tong.sbatch ../${barcode}.fastq
else
    echo "Analysis has already been run with Tong's Pipeline."
fi

