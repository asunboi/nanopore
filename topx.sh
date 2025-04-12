#!/bin/bash

awk '{print $2}' final_reads_count.txt | head -n $1 > final_readnames.txt
sed 's/^/>/' final_readnames.txt > final_readnames_grep.txt
grep -A 1 -f final_readnames_grep.txt original.fa | grep -v -- "^--$" > top$1full.fa
grep -A 1 -f final_readnames_grep.txt RHS.fa | grep -v -- "^--$" > top$1RHS.fa
