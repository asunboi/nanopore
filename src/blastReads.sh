#!/bin/bash

query=$1
iteration=$2
strand=$3

blastn -import_search_strategy ~/nanopore/search_strategy.asn -query $query -subject ~/nanopore/modRepair.fasta -strand $strand \
-outfmt "6 qseqid qstart qend" 2>/dev/null | sort -k1,1 -k2,2n -u > ${iteration}_filter.txt

cut -f1 ${iteration}_filter.txt | sort -u > ${iteration}_readnames.txt

seqtk subseq -t $query ${iteration}_readnames.txt > ${iteration}_reads.tsv
