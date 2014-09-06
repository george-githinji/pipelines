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

identities=(98 96 94 92 90 88)

seqfilefullpath=$( cd $(dirname $seqfile); pwd)/$(basename $seqfile)

#create the directories to hold output for each cluster identity
for id in ${identities[@]}; do
  mkdir "${id}"
  fr_id=$(bc<<<"scale=2; $id/100")
  echo $fr_id
  cd-hit-est -i $seqfile -o $id/cluster.${id} -c $fr_id

  cd "${id}"

  #Write a file containing the cluster-id and respective list of cluster members
  echo "Parsing CDHIT clusters"

  ruby  ~/Softwares/parse-cdhit.rb cluster.$id.clstr member >seqs_per_cluster.txt
  ruby ~/Softwares/parse-cdhit.rb cluster.$id.clstr size >seq_numbers_per_cluster.txt

  echo "Sorting cluster numbers"
  sort -t$',' -k2 -nr seq_numbers_per_cluster.txt >sorted.seq_numbers_per_cluster.txt
  awk -F, '{if($2 > 1) print $1}'< sorted.seq_numbers_per_cluster.txt | sort -nk1 >most_clusters.txt

  #create directories for each cluster in most_cluster.txt
  while read clustr; do
    mkdir cluster_$clustr
    cd cluster_$clustr
      
      echo "writing fasta file for cluster_$clustr"
      #write the fasta members for clustr $name
      pcregrep "^$clustr:" ../seqs_per_cluster.txt | awk -F: '{print $2}' | tr , '\n' |
      awk '{gsub("_","\\_",$0);$0="(?s)^>"$0".*?(?=\\n(\\z|>))"}1' | pcregrep -oM -f - $seqfilefullpath >cluster_$clustr.fasta
      
      echo "Translating"
      #translate the DNA file
      ruby ~/Softwares/translate.rb -i cluster_$clustr.fasta -f 1 -o cluster_$clustr.aa.fasta
     
      echo "Aligning sequences"
      #align both DNA and AA sequences
      muscle -in cluster_$clustr.fasta -out cluster_$clustr.aln
      muscle -in cluster_$clustr.aa.fasta -out cluster_$clustr.aa.aln
      
      echo "finding mutations from the alignments"
      #find the mutations 
      ruby ~/Softwares/mutation.rb cluster_$clustr.aln >cluster_$clustr.mutations.txt

      echo "listing mutation positions"
      #list the positions for each mutation
      awk -F, '{print $1}' cluster_$clustr.mutations.txt | sort -nrk 1 |
      uniq  >cluster_$clustr.mutation.positions.txt

      echo "Sanitizing"
       sed 's/A-G/G-A/g; s/T-C/C-T/g; s/A-T/T-A/g; s/A-C/C-A/g; s/C-G/G-C/g; s/T-G/G-T/g; s/--T/T--/g; s/--G/G--/g; s/--C/C--/g; s/--A/A--/g;' <cluster_$clustr.mutations.txt | 
       awk -F, '{print $1,$2,$3,$4}' | sort -nrk 1 | uniq >cluster_$clustr.mutations.edited.txt

       echo "Listing mutations per codon position"
       awk '{print $2,$4}' <cluster_$clustr.mutations.edited.txt | sort | uniq -c >cluster_$clustr.mutations_percodon_position.txt
       
       echo "printing dN-dS changes"
       ~/Softwares/selection.rb -i cluster_$clustr.aln -S | awk '{print $4}' | sort | uniq -c >cluster_$clustr.dN_dS.txt
    cd ..
  done <most_clusters.txt

  echo "Combining all the mutations for $id identity threshold"
  find . -type f -name *.mutations.edited.txt | xargs cat | awk '{print $2,$4}'|sort | uniq -c >all.mutations.txt

  cd ..
done
