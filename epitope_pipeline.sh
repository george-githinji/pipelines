#!/bin/sh
## This script should be improved with GNU parallel! :) 

prediction_method="NetMHCII-pan-3"

fasta_file=$1
allele_file=$2
peptide_len=$3
#output_file=$4

echo "Predicting class II T-cell epitopes" #$file_no" #>> $LOG_FILE

function wait_run_in_parallel()
{
  local number_to_run_concurrently=$1
  if [ `jobs -np | wc -l` -gt $number_to_run_concurrently ]; then
    wait `jobs -np | head -1` # wait for the oldest one to finish
  fi
}

#run the prediction and sanitize the output to csv file
function predict(){
  local allele=$1
  echo "running $allele"
  netMHCII-pan-3 -f $fasta_file -a $allele -length $peptide_len | 
  sed '/\#/d' | sed '/---/d' |
  sed '/Pos/d' | sed '/Number/d' | 
  sed 's/^ *//' | sed 's/[ \t]*$//' |
  sed 's/<=WB//' | sed 's/<=SB//' |
  tr -s '[:blank:]' ',' |
  sed '/^$/d'| sed 's/,$//'
}

# function to run MHCII program for each alleles and send the jobs to the background
run_prediction () {
  #Read the allele file
  while read mhc_allele
  do
    predict $mhc_allele >>"$mhc_allele.pred.temp" &
    wait_run_in_parallel 10 #run 10 jobs at a time and wait to execute.
  done <$allele_file
}

#execute the jobs
run_prediction
#wait for all background jobs
wait

#concatenate the files
cat *.pred.temp >>all.hla.preds

#cleanup!
rm *.pred.temp
