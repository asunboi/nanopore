#!/bin/bash

mkdir cleaning

dir_name=`pwd`
barcode=${dir_name##*/}

cat *.gz > cleaning/${barcode}.fastq.gz
cd cleaning

gunzip ${barcode}.fastq.gz

## forward reads
cutadapt -g XGATAATGATCATAATCAGCCATACC -e 0.2 --overlap 10 --discard-untrimmed -o forward.fa ${barcode}.fastq

## reverse reads
cutadapt -a GGTATGGCTGATTATGATAGTTATCX -e 0.2 --overlap 10 --discard-untrimmed -o reverse.fa ${barcode}.fastq

## reads that somehow have both on same strand (unlikely, but we cannot differentiate orientation so we will have to throw away)
cutadapt -g XGATAATGATCATAATCAGCCATACC...GGTATGGCTGATTATGATAGTTATCX -e 0.2 --overlap 10 --discard-untrimmed -o dupes.fa ${barcode}.fastq

grep ">" dupes.fa | sed 's/^.//'| cut -d " " -f 1 > dupes.txt

faSomeRecords -exclude forward.fa dupes.txt forwardND.fa
faSomeRecords -exclude reverse.fa dupes.txt reverseND.fa

seqtk seq -r reverseND.fa > flipReverseND.fa

cat forward.fa > tmp.fa
cat flipReverseND.fa >> tmp.fa

#removes the leftmost instance of the reverse adaptor in the correct configuration
cutadapt -a AGATCGGAAGAGCACACGTCTG -e 0.2 --overlap 10 -o noprimer.fa tmp.fa

### blasting and removing forward reporter portions ###

blastReads.sh noprimer.fa 1 'plus'
python ~/scripts/nanopore/cutReads.py 1_reads.tsv 1_filter.txt 1_filter.fa

iteration=1

while :
do
        nextIteration=$(($iteration + 1))
        blastReads.sh ${iteration}_filter.fa $nextIteration 'plus'
        if [[ ! -s ${nextIteration}_filter.txt ]];
        then
                #statements
                echo "no more forward reporter regions"
                appendReads.sh ${iteration}_filter.fa ${nextIteration}_readnames.txt
                break
        else
                appendReads.sh ${iteration}_filter.fa ${nextIteration}_readnames.txt
                python ~/scripts/nanopore/cutReads.py ${nextIteration}_reads.tsv ${nextIteration}_filter.txt ${nextIteration}_filter.fa
                iteration=$(($iteration + 1))
        fi
done

### blasting and removing reverse reporter portions ###

while :
do
        nextIteration=$(($iteration + 1))
        blastReads.sh ${iteration}_filter.fa $nextIteration 'minus'
        if [[ ! -s ${nextIteration}_filter.txt ]];
        then
                #statements
                echo "no more reverse reporter regions"
                appendReads.sh ${iteration}_filter.fa ${nextIteration}_readnames.txt
                break
        else
                appendReads.sh ${iteration}_filter.fa ${nextIteration}_readnames.txt
                python ~/scripts/nanopore/cutReads.py ${nextIteration}_reads.tsv ${nextIteration}_filter.txt ${nextIteration}_filter.fa
                iteration=$(($iteration + 1))
        fi
done

