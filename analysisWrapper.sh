#!/bin/bash

mkdir $1
mkdir $1/cleaning
mkdir $1/clustering
mkdir $1/tong

cd $1
sbatch ~/scripts/nanopore/analysis.sbatch $1

