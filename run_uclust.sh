#!/bin/bash

#Check to see if at least one argument was specified
if [ $# -lt 1 ] ; then
  echo "You must specify the directory -d and identity -i arguments"
  exit 1
fi

#Process the arguments
while getopts d:i: opt
do
  case "$opt" in
    d) dir=$OPTARG;;
    i) identity=${OPTARG};;
esac
done
 echo $dir
 echo $identity

for file in $dir/*.fasta
  do
  filename=$(basename $file)
  filename="${filename%.*}"
  
  #sort the sequence
  uclust --sort $file --output $filename.sorted.fasta

  #cluster at given identity
  uclust --input $filename.sorted.fasta --uc $filename.uc --id ${identity}

  #write a fasta output of all the seed sequences as a nr database
  uclust --input $filename.sorted.fasta --uc2fasta $filename.uc --types S --output $filename.$identity.fasta
done

mkdir -p sorted_fasta
mv *.sorted.fasta sorted_fasta

mkdir -p cluster_files 
mv *.uc cluster_files

#sort the cluster file
sort -n -k2 -k4 cluster_files/*.uc >cluster_files/sorted.nr.uc

#from the sorted cluster file write a new file containing the frequency of members for each cluster
awk '{print $2}' sorted.nr.uc | sort -n|uniq -c | awk '{if($1>2) print}' | sort -nk1 >clusters.names.sizes.txt
