#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24gb
#SBATCH --time=24:00:00
    
module load blast/2.13.0+

grep "^>" RHS.fa > RHSnames.txt
grep -A 1 -f RHSnames.txt LHS.fa | grep -v -- "^--$" > LHS_filtered.fa

# Make a blast database, then blast it against itself (All vs All)
makeblastdb -in LHS_filtered.fa -dbtype nucl -parse_seqids
blastn -import_search_strategy ~/nanopore/short_search_strategy.asn -query LHS_filtered.fa -db LHS_filtered.fa -num_threads 8 \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
-out LHS.blast

echo blasting generate `wc -l LHS.blast | awk '{print $1}'` alignment in total

# Co-coverage >= to 80% and Percent Identity is bigger than 80%, the read passes the filter. 

awk '{
    if ($3 > 80 && ($8 >= ($13 - 50) && $10 >= ($14 - 50))) {
        print;
    }
}' LHS.blast > filtered_LHS.blast

echo `wc -l filtered_LHS.blast | awk '{print $1}'` LHS alignment left after filtering

# Make a blast database, then blast it against itself (All vs All)
makeblastdb -in RHS.fa -dbtype nucl -parse_seqids
blastn -import_search_strategy ~/nanopore/short_search_strategy.asn -query RHS.fa -db RHS.fa -num_threads 8 \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
-out RHS.blast

echo blasting generate `wc -l RHS.blast | awk '{print $1}'` alignment in total

# Co-coverage >= to 80% and Percent Identity is bigger than 85%, the read passes the filter. 

awk '{
    if ($3 > 80 && ($7 <= 50 && $9 <= 50)) {
        print;
    }
}' RHS.blast > filtered_RHS.blast

echo `wc -l filtered_RHS.blast | awk '{print $1}'` RHS alignment left after filtering

sort -k1,1 -k14,14nr RHS_matching_alignments.txt > sorted_RHS_matching_alignments.txt
awk 'NR==1 {prev=$1; print; next} $1 != prev {print} {prev=$1}' sorted_RHS_matching_alignments.txt > long_RHS_alignments.txt
sort -k2,2 long_RHS_alignments.txt > sorted_long_RHS_alignments.txt
awk '{print $2}' sorted_long_RHS_alignments.txt | uniq -c > sorted_long_RHS_alignments_count.txt
sort -r -n k1,1 sorted_long_RHS_alignments_count.txt > final_reads_count.txt

