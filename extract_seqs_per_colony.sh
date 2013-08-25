#!/bin/bash

#PARSED_OPTIONS=$(getopt -n "$0" -0 h sqc: --long "seqfile,qualfile,cloneslist:" -- "$@")
##bad arfgs

#if [ $? -ne 0];
#then
  #exit 1
#fi

#eval set -- "$PARSED_OPTIONS"

#while true;:w

#do
  #case "$1" in 
    #-h|--help)
      #shift;;

      #-s|--seq)

  #esac
#done

while read clone
do
    cat seqs/BGI_reads.fasta | perl -e '$ok=0; while(<STDIN>){if(/^>/&&/'$clone'/){$ok=1;print} elsif(/^>/){$ok=0}elsif($ok){print}}' >data_per_colony/$clone.fasta
    cat seqs/BGI_reads.qual  | perl -e '$ok=0; while(<STDIN>){if(/^>/&&/'$clone'/){$ok=1;print} elsif(/^>/){$ok=0}elsif($ok){print}}' >data_per_colony/$clone.fasta.qual
done <clones_list.txt
