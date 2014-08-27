#!/bin/sh

prediction_method="NetMHCII-pan-3"

fasta_file=$1
allele_file=$2
output_file=$3

echo "Running compute peptides for file" #$file_no" #>> $LOG_FILE


# function to run MHCII program with different alleles
run_prediction () {

#Read the allele file
while read mhc_allele
do

  echo "BINDING_TEST $mhc_allele " >>$output_file

  NetMHCIIpan3 -f $fasta_file $mhc_allele >>$output_file || exit e

done < $allele_file

}
