#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24gb
#SBATCH --time=24:00:00
    
module load blast/2.13.0+

# Make a blast database, then blast it against itself (All vs All)
makeblastdb -in trimmed.fasta -dbtype nucl -parse_seqids
blastn -import_search_strategy ~/nanopore/short_search_strategy.asn -query trimmed.fasta -db trimmed.fasta -num_threads 8 \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
-out seq.blast

echo blasting generate `wc -l seq.blast | awk '{print $1}'` alignment in total

# If Query length < Sequence length, min = Query length. Else, min = sequence length. 
# If alignment length is >= to 90% of the minimum length and Percent Identity is bigger than 85%, the read passes the filter. 
awk '{if($13 < $14) {min=$13} else {min=$14}}{if(($4 >= (min*0.9)) && ($3 >= 85)){print $0}}' seq.blast > seq_filtered.blast

echo `wc -l seq_filtered.blast | awk '{print $1}'` alignment left after filtering

awk '{print $1}' seq_filtered.blast | uniq -c > query.count
start=1
while read line; do echo $line $start; id=($line); let start=start+${id[0]}; done < query.count > query.line.count
sort -g -r query.line.count > query.sorted.line.

######split large alignment file into smaller scale (1M lines per file) to accelerate later processing speed
echo split large alignment file into smaller scale 1M lines per file
n=1
line_num=1
cutoff=1
lines=`wc -l query.count | awk '{print $1}'`
for i in `seq 1 $lines`
do id=(`awk 'NR == '$i'{print; exit}' query.count`)
echo ${id[0]} ${id[1]} $line_num >> ${n}.temp.line.count
let line_num=line_num+${id[0]}
if [ $line_num -gt 1000000 ]
then let upper=$line_num+cutoff-1
awk '{if (NR >= '$cutoff'){if (NR < '$upper') {print $0} else {exit}}}' seq_filtered.blast > ${n}.temp.blast
cutoff=$upper
let n=n+1
line_num=1
fi
if [ $i -eq $lines ]
then awk 'NR >= '$cutoff'' seq_filtered.blast > ${n}.temp.blast
fi
done

######take all query which has more than 20 alignments, and take the longest read in each pool as representative to make a pool only containing long reads
    echo filter count file to exclude reads with less than 50 counts and take the long representative read as first round of de-scaling
temp_num=`ls *.temp.blast | wc -l`
for i in `seq 1 $temp_num`
do while read line
do id=($line)

## FILTERING STEP
if [ ${id[0]} -lt 50 ]
then echo ${id[0]} ${id[1]} >> fail.list
continue
fi
awk '{if (NR >= '${id[2]}'){if (NR < ('${id[2]}'+'${id[0]}')) {print $0} else {exit}}}' ${i}.temp.blast > temp.pool
sort -g -r -k 14 temp.pool > temp.sorted.pool
head -n 1 temp.sorted.pool >> long.pool
done < ${i}.temp.line.count
done
awk '{print $2}' long.pool | sort | uniq -c | sort -g -r > long.count

######use the read id from long count file to recover sequence from fasta file and do the second alignemnt and filtering
echo second round of de-scaling
for seq in `awk '{print $2}' long.count | cat`
do grep -A1 $seq trimmed.fasta
done > second.fasta
makeblastdb -in second.fasta -dbtype nucl -parse_seqids
blastn -import_search_strategy ~/nanopore/short_search_strategy.asn -query second.fasta -db second.fasta -num_threads 8 \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
-out second.blast
echo blasting generate `wc -l second.blast | awk '{print $1}'` alignment in total
awk '{if($13 < $14) {min=$13} else {min=$14}}{if(($4 >= (min*0.9)) && ($3 >= 85)){print $0}}' second.blast > second.filtered.blast
echo `wc -l second.filtered.blast | awk '{print $1}'` alignment left after filtering
awk '{print $1}' second.filtered.blast | uniq -c > second.query.count
start=1
while read line; do echo $line $start; id=($line); let start=start+${id[0]}; done < second.query.count > second.query.line.count
######take long pool from the second blast pool and count, all remaining reads together can represent all original sequence in the beginning
echo filter the second count file to exclude reads with less than 10 counts
while read line
do id=($line)
if [ ${id[0]} -lt 10 ]
then count=`grep ${id[1]} long.count | awk '{print $1}'`
if [ $count -lt 50 ]
then echo pool ${id[0]} smaller than 50
echo $line $count >> second.fail.list
continue
fi
fi
awk '{if (NR >= '"${id[2]}"'){if (NR < ('"${id[2]}"'+'"${id[0]}"')) {print $0} else {exit}}}' second.filtered.blast > temp.pool
sort -g -r -k 14 temp.pool > temp.sorted.pool
head -n 1 temp.sorted.pool >> second.long.pool
done < second.query.line.count
awk '{print $2}' second.long.pool | sort | uniq -c | sort -g -r > second.long.count


######grouping the reads in second count file into final pools
echo grouping final long reads into individual pools
line_num=`wc -l second.long.count| awk '{print $1}'`
for i in `seq 1 $line_num`
do id=(`awk 'NR == '$i' {print; exit}' second.long.count`)
if [ $i == 1 ]
then let final_pool_count=1
echo ${id[1]} > ${final_pool_count}.final.pool
continue
fi
for j in `seq 1 $final_pool_count`
do while read line
do grepF=`grep ${id[1]}$'\t'$line second.filtered.blast`
grepR=`grep $line$'\t'${id[1]} second.filtered.blast`
if [[ $grepF ]] || [[ $grepR ]]
then echo ${id[1]} can align to $line
echo ${id[1]} >> ${j}.final.pool
break 2
fi
done < ${j}.final.pool
if [ $j == $final_pool_count ]
then echo ${id[1]} represent a new pool
let final_pool_count=final_pool_count+1
echo ${id[1]} > ${final_pool_count}.final.pool
fi
done
done


######put all sequence together into seperate file
######sometimes sequence in long count may be filtered out due to low alignment number in the very first step
pool_num=`ls *final.pool | wc -l`
for i in `seq 1 $pool_num`
do while read line
do grep $'\t'$line$'\t' second.long.pool
done < ${i}.final.pool > ${i}.final.sequence
awk '{print $1}' ${i}.final.sequence | sort | uniq > ${i}.temp.final.pool
awk '{print $2}' ${i}.final.sequence | uniq >> ${i}.temp.final.pool
sort ${i}.temp.final.pool | uniq > ${i}.final.pool
while read line
do grep $'\t'$line$'\t' long.pool
done < ${i}.final.pool > ${i}.final.sequence
awk '{print $1}' ${i}.final.sequence | sort | uniq > ${i}.temp.final.pool
awk '{print $2}' ${i}.final.sequence | uniq >> ${i}.temp.final.pool
sort ${i}.temp.final.pool | uniq > ${i}.final.pool
rm ${i}.temp.final.pool
done


######calculate how many reads in each pool according to the first long pool and store alignment information in corresponding file
echo counting exact read number from the first long count file
for i in `seq 1 $pool_num`
do while read line
do grep $line long.count
done < ${i}.final.pool | sort -g -r > ${i}.final.pool.count
while read line
do id=($line)
grep $'\t'${id[1]}$'\t' -m 1 long.pool | awk '{print $2, $14}'
grep $'\t'${id[1]}$'\t' long.pool | awk '{print $1, $13}'
done < ${i}.final.pool.count | sort | uniq | sort -g -r -k 2 > ${i}.final.sequence
done


######count number for each pool and representative sequence will be stored in final.sorted.count
dir_name=`pwd`
prefix=${dir_name##*/}
echo count number for each pool and representative sequence will be stored in ${prefix}.sorted.count

for i in `seq 1 $pool_num`
do count1=`awk '{sum+=$1}END{print sum}' ${i}.final.pool.count`
count2=`wc -l ${i}.final.sequence | awk '{print $1}'`
######two different ways to calculate read count, there should be no big difference between them.
######some sequences may present in long pool but itself is filtered out due to less than 50 alignments.
if [[ $count1 -lt 50 ]] && [[ $count2 -lt 50 ]]
then continue
fi
rep=(`head -n 1 ${i}.final.sequence`)
seq=(`grep -A1 ${rep[0]} original.fa | tail -n 1`)
echo ">"${rep[0]} $count1 $count2 ${i}.final.pool ${rep[1]}
echo $seq
done > ${prefix}.rep.fasta
