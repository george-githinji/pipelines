#!/bin/bash

#This script should be run from the  grouped_by_clones folder.
#i.e all the sequences for each isolate's clones are located
#please run the copy_per_colony script first to generate this list of sequences

#get the list of isolates from ls the current folder
echo "generating a list of isolate names"

sed '/^$/d' isolate_names.txt >isolate_names
find . -type f -name '*.fasta' | grep -Eo '[0-9]{1,}[A-Z]{1}' | sort | uniq >isolate_names.txt

while read isolate
do
  mkdir -p ../$isolate

  echo "copying data for $isolate"
  cp $isolate* ../$isolate

  cd ../$isolate

  ls *.fasta | perl -lne 'print /.{3}-?(.+)\..+[^\>]/' >$isolate\_clones.txt

  #remove empty lines
  sed '/^$/d' $isolate\_clones.txt | uniq >temp.txt
  mv temp.txt $isolate\_clones.txt

  while read clone
  do
    echo "creating clone $clone directory"

    mkdir -p $isolate\_$clone

    echo "copying $isolate-$clone.fasta and $isolate-$clone.fasta.qual files to directory $clone"

    mv $isolate$clone\.* $isolate\_$clone

    cd $isolate\_$clone

    echo "running cap3 for $clone read pairs"
    isolatename="$isolate$clone" 
    ~/cap3/cap3 $isolatename\.fasta $isolatename\.fasta.qual >$isolate\_$clone\.cap.txt

    #append the clone names to the contigs rename the contigs 
    sed "s/^>/>$isolate\-$clone\_/" $isolatename\.fasta.cap.contigs >$isolate\_$clone\.renamed.contigs.fasta
    sed "s/^>/>$isolate\-$clone\_/" $isolatename\.fasta.cap.contigs.qual >$isolate\_$clone\.renamed.contigs.qual

    #create a single line fasta and qual for the assembles and renamed contigs
    awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $isolate\_$clone\.renamed.contigs.fasta >$isolate\_$clone\.SL.renamed.contigs.fasta 
    awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $isolate\_$clone\.renamed.contigs.qual >$isolate\_$clone\.SL.renamed.contigs.qual 

    sed '/^$/d' $isolate\_$clone\.SL.renamed.contigs.qual >temp.txt
    mv temp.txt $isolate\_$clone\.SL.renamed.contigs.qual

    sed '/^$/d' $isolate\_$clone\.SL.renamed.contigs.fasta >temp.txt
    mv temp.txt $isolate\_$clone\.SL.renamed.contigs.fasta

    #create a fastq file from the assembled contig and the corresponding quality file
    echo "creating a fastq file"
    #perl ~/Softwares/makefastq.pl $isolate\_$clone\.SL.renamed.contigs.fasta $isolate\_$clone\.SL.renamed.contigs.qual >$isolate\_$clone\.SL.renamed.contigs.fastq
    python ~/Softwares/to_fastq.py $isolate\_$clone\.SL.renamed.contigs.fasta $isolate\_$clone\.SL.renamed.contigs.qual $isolate\_$clone\.SL.renamed.contigs.fastq


    echo "calling dbla tags for $clone contig"
    ruby ~/Code/Ruby/bioruby-dbla-finder/lib/bio-dbla-finder.rb $isolate\_$clone\.SL.renamed.contigs.fastq >$isolate\_$clone\.SL.tag.renamed.contigs.fastq
    #dbla_finder -i $isolate\_$clone\.SL.renamed.contigs.fastq >$isolate\_$clone\.SL.tag.renamed.contigs.fasta

    cd ..

  done <$isolate\_clones.txt

  cd ../data_per_colony

done <isolate_names

cd ..

#echo "concatenate the fasta files"

#find . -type f -name '*.SL.tag.renamed.contigs.fastq' -print0 | xargs -0 cat >tags.fastq

#print tags.fasta
#ruby ~/Softwares/read_fastq.rb tags.fastq seq >tags.fasta

#count the total number of potential tags
#countseqs tags.fasta

#sometimes the dblfinder detects multitple tag like regions in longer contigs.
#These need to be corrected/removed manually
#grep ">" tags.fasta | sort | uniq -c | sort -nk 1 | uniq | awk '{if($1>1) print}' >duplicated_tags.txt
