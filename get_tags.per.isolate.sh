#!/bin/bash

#read the isolates names

while read isolate
do
  echo $isolate
  cd ../$isolate
find . -type f -name '$isolate*.tag.renamed.contigs.fastq' -print0 | xargs -0 cat #>$isolate.tags.fastq
#ruby ~/Softwares/read_fastq.rb $isolate.tags.fastq >$isolate.tags.fasta
  cd ..
done<$1
