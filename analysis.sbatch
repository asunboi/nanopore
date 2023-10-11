#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24gb
#SBATCH --time=24:00:00

module load python

mkdir $1
mkdir $1/cleaning
mkdir $1/clustering
mkdir $1/tong

### check whether the compiled fasta has already been made - if not, make it.
dir_name=`pwd`
barcode=${dir_name##*/}

if ! test -f ${barcode}.fa; then
  cat *.gz > ${barcode}.fastq.gz
  gunzip ${barcode}.fastq.gz
  seqtk seq -a ${barcode}.fastq > ${barcode}.fa
fi

cd $1/tong 

sbatch ~/scripts/nanopore/tong.sbatch ../../${barcode}.fastq

cd ..

cd cleaning

## forward reads
cutadapt -g XGATAATGATCATAATCAGCCATACC -e 0.2 --overlap 10 --discard-untrimmed -o forward.fa ../../${barcode}.fastq

## reverse reads
cutadapt -a GGTATGGCTGATTATGATAGTTATCX -e 0.2 --overlap 10 --discard-untrimmed -o reverse.fa ../../${barcode}.fastq

## reads that somehow have both on same strand (unlikely, but we cannot differentiate orientation so we will have to throw away)
cutadapt -g XGATAATGATCATAATCAGCCATACC...GGTATGGCTGATTATGATAGTTATCX -e 0.2 --overlap 10 --discard-untrimmed -o dupes.fa ../../${barcode}.fastq

grep ">" dupes.fa | sed 's/^.//'| cut -d " " -f 1 > dupes.txt

faSomeRecords -exclude forward.fa dupes.txt forwardND.fa
faSomeRecords -exclude reverse.fa dupes.txt reverseND.fa

seqtk seq -r reverseND.fa > flipReverseND.fa

cat forward.fa > tmp.fa
cat flipReverseND.fa >> tmp.fa

#removes the leftmost instance of the reverse adaptor in the correct configuration (4x)
cutadapt -a AGATCGGAAGAGCACACGTCTG -e 0.2 --overlap 10 -n 4 -m 1 -o noprimer.fa tmp.fa

### scuffed symlink, TODO
cp tmp.fa ../clustering/original.fa

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
                appendReads.sh ${iteration}_filter.fa ${nextIteration}_readnames.txt 'noForwardRep.fa'
                break
        else
                appendReads.sh ${iteration}_filter.fa ${nextIteration}_readnames.txt 'noForwardRep.fa'
                python ~/scripts/nanopore/cutReads.py ${nextIteration}_reads.tsv ${nextIteration}_filter.txt ${nextIteration}_filter.fa
                iteration=$(($iteration + 1))   
        fi
done

### blasting and removing reverse reporter portions ###

inputFile='noForwardRep.fa'

while :
do
        nextIteration=$(($iteration + 1))
        blastReads.sh $inputFile $nextIteration 'minus'
        if [[ ! -s ${nextIteration}_filter.txt ]];
        then
                #statements
                echo "no more reverse reporter regions"
                appendReads.sh $inputFile ${nextIteration}_readnames.txt 'noReverseRep.fa'
                break
        else
                appendReads.sh $inputFile ${nextIteration}_readnames.txt 'noReverseRep.fa'
                python ~/scripts/nanopore/cutReads.py ${nextIteration}_reads.tsv ${nextIteration}_filter.txt ${nextIteration}_filter.fa
                inputFile=${nextIteration}_filter.fa
                iteration=$(($iteration + 1))
        fi
done

## output is noReverseRep.fa

RepeatMasker -noint noReverseRep.fa

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < noReverseRep.fa.masked | sed 's/N//g' | tr "\t" "\n" > noSimple.fa

seqtk seq -L 50 noSimple.fa > cleaned.fa

cd ..

cp cleaning/cleaned.fa clustering/trimmed.fasta

cd clustering

sbatch ~/scripts/nanopore/cluster.sbatch