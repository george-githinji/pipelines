#!/bin/sh

prediction_method="NetMHCII-pan-3"

fasta_file=$1
allele_file=$2
peptide_len=$3
output_file=$4

echo "Predicting class II T-cell epitopes" #$file_no" #>> $LOG_FILE

# function to run MHCII program with different alleles
run_prediction () {
  #Read the allele file
  while read mhc_allele
    do
      echo "$mhc_allele"
      netMHCII-pan-3 -f $fasta_file -a $mhc_allele -length $peptide_len | 
      sed '/\#/d' | sed '/---/d' |
      sed '/Pos/d' | sed '/Number/d' | 
      sed 's/^ *//' | sed 's/[ \t]*$//' |
      sed 's/<=WB//' | sed 's/<=SB//' |
      tr -s '[:blank:]' ',' |
      sed '/^$/d'| sed 's/,$//' >>$output_file
    done <$allele_file
}

run_prediction
