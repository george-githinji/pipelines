#!/bin/bash
#provide consitent naming convection

for file in *.fasta
do
echo "processing $file"
sed 's/_Contig1//' fasta/$file | sed 's/B//' | sed 's/-/./' >$file.fasta
cd-hit-est 
done

