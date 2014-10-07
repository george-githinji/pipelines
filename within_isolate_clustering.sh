#!/bin/bash

#find the isolate file name in seqs folder
for file in ../seqs/*.fa
do
  fullpath=$(cd $(dirname $file); pwd)/$(basename $file)
  filepath="${file##*/}"
  filename="${filepath%.*}"
  mkdir -p $filename
  cd $filename
  ruby ~/Softwares/get_blocks.rb $fullpath >$filename.block_H.fasta
  clusteringpipeline $filename.block_H.fasta
  cd ..
done

mkdir block_H
mv *.nr block_H
