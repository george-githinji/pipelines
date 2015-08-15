#!/bin/bash
DIR=$1

for file in $DIR/*.fa
do
  filename=$(basename $file)
  ruby ~/Code/Ruby/blocks/lib/polymorphic_sites.rb $file >$filename.snps.txt
done
