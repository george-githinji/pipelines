#!/bin/bash

seq_file=$1

filename=$(basename $seq_file)

#create a blastdb of the sequences
makeblastdb -in $seq_file -dbtype nucl -out $filename.blastdb

#run all vs all blast against the database with 4 threads
blastn -db $filename.blastdb -query $seq_file -outfmt 6 -out all-vs-all.tsv -num_threads 4

#sequence identity > 90% and the alignment length > 300
awk '{if($1 != $2 && $3 > 88 && $4 > 300) print $1,$2,$3}' all-vs-all.tsv | awk '!seen[$1<=$2? $1 FS $2: $2 FS $1]++ {print $1,$2,$3}'  >all-vs-all.88.abc
awk '{if($1 != $2 && $3 > 90 && $4 > 300) print $1,$2,$3}' all-vs-all.tsv | awk '!seen[$1<=$2? $1 FS $2: $2 FS $1]++ {print $1,$2,$3}'  >all-vs-all.90.abc
awk '{if($1 != $2 && $3 > 92 && $4 > 300) print $1,$2,$3}' all-vs-all.tsv | awk '!seen[$1<=$2? $1 FS $2: $2 FS $1]++ {print $1,$2,$3}'  >all-vs-all.92.abc
awk '{if($1 != $2 && $3 > 94 && $4 > 300) print $1,$2,$3}' all-vs-all.tsv | awk '!seen[$1<=$2? $1 FS $2: $2 FS $1]++ {print $1,$2,$3}'  >all-vs-all.93.abc
awk '{if($1 != $2 && $3 > 96 && $4 > 300) print $1,$2,$3}' all-vs-all.tsv | awk '!seen[$1<=$2? $1 FS $2: $2 FS $1]++ {print $1,$2,$3}'  >all-vs-all.96.abc
awk '{if($1 != $2 && $3 > 98 && $4 > 300) print $1,$2,$3}' all-vs-all.tsv | awk '!seen[$1<=$2? $1 FS $2: $2 FS $1]++ {print $1,$2,$3}'  >all-vs-all.98.abc


#Usearch clustering
identities=(0.98 0.96 0.94 0.92 0.90 0.88)

for id in ${identities[@]}
do
  vsearch -cluster_smallmem $seq_file -id $id -centroids $filename.$id.nr.fasta -uc $filename.$id.uc -blast6out $filename.$id.m6
done

#run mcl with different inflation parameters
#for i in {1..10}
#   do mcl $all-vs-all.abc --abc -I $i -te 16 -o ${base}_$i.out
#done
