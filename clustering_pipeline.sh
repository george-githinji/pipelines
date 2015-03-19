#!/bin/bash
#######################################################
#This bash script takes a sequence file and clusters the sequences at predifined identities (98 - 88)
#For each cluster in each identity, we align the sequences and search for changes, positions of the changes and the type of changes. 
#We tally all the changes in each identity

#Dependencies
#Bio-cd-hit-report
#Ruby
#AWK
#sed
#translate.rb
#mutation.rb

#Author: George Githinji
#Email: ggithinji@kemri-wellcome.org
#MIT Licence
#######################################################

seqfile=$1

[ $# -eq 0 ] && { echo "Usage: $0 nucleotide_file"; exit 1;}

identities=(98 96 94 92 90 88)

seqfilefullpath=$( cd $(dirname $seqfile); pwd)/$(basename $seqfile)

#sort the input sequence with usearch sortbylength command
vsearch -sortbylength $seqfile --output $seqfile.sorted.fasta

sortedfilefullpath=$( cd $(dirname $seqfile); pwd)/$(basename $seqfile.sorted.fasta)

#create the directories to hold output for each cluster identity
for id in ${identities[@]}; do
  mkdir -p "${id}"
  fr_id=$(bc<<<"scale=2; $id/100")
  echo $fr_id

  #use cd-hit to cluster
  #cd-hit-est -i $seqfile -o $id/cluster.${id} -c $fr_id #retire cd-hit

  #use the uclust algorithm to perform the clustering
  #uclust --input $sortedfilefullpath --uc $id/cluster.${id}.uc --id ${fr_id} --optimal

  #use vsearch1.1.1
  vsearch -cluster_smallmem $sortedfilefullpath -id ${fr_id} -centroids $id/cluster.${id}.nr.fasta -uc $id/cluster.${id}.uc

  #convert the output to cd-hit format
  uclust --uc2clstr $id/cluster.${id}.uc --output $id/cluster.${id}.clstr

  cd "${id}"

  #Write a file containing the cluster-id and respective list of cluster members
  #echo "Parsing clusters"

  ruby ~/Softwares/parse-cdhit.rb cluster.$id.clstr member >seqs_per_cluster.txt
  ruby ~/Softwares/parse-cdhit.rb cluster.$id.clstr size >seq_numbers_per_cluster.txt

  echo "Sorting cluster numbers"
  sort -t$',' -k2 -nr seq_numbers_per_cluster.txt >sorted.seq_numbers_per_cluster.txt
  awk -F, '{if($2 > 1) print $1}'< sorted.seq_numbers_per_cluster.txt >most_clusters.txt

  #create directories for each cluster in most_cluster.txt
  while read clustr; do
    echo "processing $clustr"

    mkdir cluster_$clustr
    cd cluster_$clustr

      #write the fasta members for clustr $name
      grep "^${clustr}:" ../seqs_per_cluster.txt | awk -F: '{print $2}' | tr , '\n' >members.txt
      fastagrep -F -X -f members.txt $sortedfilefullpath >cluster_$clustr.fasta

      #echo "Translating"
      #ruby ~/Softwares/translate.rb -i cluster_$clustr.fasta -f 1 -o cluster_$clustr.aa.fasta

      #echo "Aligning sequences"
      #muscle -in cluster_$clustr.fasta -out cluster_$clustr.aln #-physout cluster_$clustr.phylip
      prank -d=cluster_$clustr.fasta -o=cluster_$clustr -codon -F -showtree -quiet

      #replace mutation.rb with polymorphic_sites.rb
      #echo "Searching for polymorphic sites"
      #ruby ~/Code/Ruby/blocks/lib/polymorphic_sites.rb cluster_$clustr.aln >cluster_$clustr.snps.muscle.txt
      ruby ~/Code/Ruby/blocks/lib/polymorphic_sites.rb cluster_$clustr.best.fas >cluster_$clustr.snps.prank.txt

      #echo "searching dN and dS changes"
      #~/Softwares/selection.rb -i cluster_$clustr.aln -S | awk '{print $4}' | sort | uniq -c | sed -e 's/^[ \t]*//' >cluster_$clustr.dN_dS.txt
      ~/Softwares/selection.rb -i cluster_$clustr.best.fas -S | awk '{print $4}' | sort | uniq -c | sed -e 's/^[ \t]*//' >cluster_$clustr.dN_dS_prank.txt
    cd ..
  done <most_clusters.txt

  echo "collect snps in $id identity threshold"
  #find . -type f -name *.snps.muscle.txt | xargs cat | sed -e 's/^[ \t]*//' >all.snps.muscle.txt
  find . -type f -name *.snps.prank.txt | xargs cat | sed -e 's/^[ \t]*//' >all.snps.prank.txt

  cd ..
done
