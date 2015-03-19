#!/bin/bash

clusterfile=$1

clusterFileFullPath=$( cd $(dirname $clusterfile); pwd)/$(basename $clusterfile)

directory=$(dirname ${clusterFileFullPath})

 while read cluster_no; do
 
   #create a neighbour joining tree
   #muscle -maketree -in $directory/cluster_$cluster_no/cluster_$cluster_no.fasta -out cluster_$cluster_no.tree -cluster neighborjoining

   #awk '{printf("%s",$0)}'  <cluster_$cluster_no.tree >cluster_$cluster_no_2.tree

   #overwrite the original tree file
   #mv cluster_$cluster_no_2.tree cluster_$cluster_no.tree

   #get a random sequence from this cluster
   #This solution is based on reservoir sampling. It is optimal in terms of both time and space complexity and works with an arbitrarily large input file.
   #https://www.biostars.org/p/18831/
   bioawk -c fastx -v k=1 '{y=x++<k?x-1:int(rand()*x);if(y<k)a[y]=">"$name"\n"$seq}END{for(z in a)print a[z]}' $directory/cluster_$cluster_no/cluster_$cluster_no.fasta >cluster_$cluster_no.fa

   Rscript ~/Code/R/phylosim_simulation.R cluster_$cluster_no.fa $directory/cluster_$cluster_no/cluster_$cluster_no.best.dnd cluster_$cluster_no.s.fa

  ruby ~/Code/Ruby/blocks/lib/polymorphic_sites.rb cluster_$cluster_no.s.fa >cluster_$cluster_no.snps.txt

 done <$clusterfile

  #find . -type f -name *.snps.txt | xargs cat | sed -e 's/^[ \t]*//' >all.snps

