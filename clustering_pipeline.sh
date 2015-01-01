
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
usearch -sortbylength $seqfile --output $seqfile.sorted.fasta

sortedfilefullpath=$( cd $(dirname $seqfile); pwd)/$(basename $seqfile.sorted.fasta)

#create the directories to hold output for each cluster identity
for id in ${identities[@]}; do
  mkdir "${id}"
  fr_id=$(bc<<<"scale=2; $id/100")
  echo $fr_id

  #use cd-hit to cluster
  #cd-hit-est -i $seqfile -o $id/cluster.${id} -c $fr_id #retire cd-hit

  #use the uclust algorithm to perform the clustering
  #uclust --input $sortedfilefullpath --uc $id/cluster.${id}.uc --id ${fr_id} --optimal

  #use usearch7
  usearch -cluster_smallmem $sortedfilefullpath -id ${fr_id} -centroids $id/cluster.${id}.nr.fasta -uc $id/cluster.${id}.uc

  #convert the output to cd-hit format
  uclust --uc2clstr $id/cluster.${id}.uc --output $id/cluster.${id}.clstr

  cd "${id}"

  #Write a file containing the cluster-id and respective list of cluster members
  echo "Parsing clusters"

  ruby ~/Softwares/parse-cdhit.rb cluster.$id.clstr member >seqs_per_cluster.txt
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
      grep "^${clustr}:" ../seqs_per_cluster.txt | awk -F: '{print $2}' | tr , '\n' >members.txt
      fastagrep -F -X -f members.txt $sortedfilefullpath >cluster_$clustr.fasta

      echo "Translating"
      #translate the DNA file
      ruby ~/Softwares/translate.rb -i cluster_$clustr.fasta -f 1 -o cluster_$clustr.aa.fasta

      echo "Aligning sequences"
      #align both DNA and AA sequences
      muscle -in cluster_$clustr.fasta -out cluster_$clustr.aln
      muscle -in cluster_$clustr.aa.fasta -out cluster_$clustr.aa.aln

     #echo "finding mutations from the alignments"

      #find the mutations 
      #this procedure takes a long time to complete.
      #ruby ~/Softwares/mutation.rb cluster_$clustr.aln >cluster_$clustr.mutations.txt

      #replace mutation.rb with polymorphic_sites.rb
      echo "Searching for polymorphic sites"
      ruby ~/Code/Ruby/blocks/lib/polymorphic_sites.rb cluster_$clustr.aln >cluster_$clustr.snps.txt

      #create a file for region coordinates
      ruby ~/Code/Ruby/blocks/lib/print_hv_positions.rb cluster_$clustr.aln >regions_positions.txt

      #split alignment to blocks and print start and stop positions
      bash ~/Code/Bash/splice_alignment.sh regions_positions.txt cluster_$clustr.aln

      #find snps in each region
      ruby ~/Code/Ruby/blocks/lib/polymorphic_sites.rb region_1.aln >snps_region1.txt
      ruby ~/Code/Ruby/blocks/lib/polymorphic_sites.rb region_2.aln >snps_region2.txt
      ruby ~/Code/Ruby/blocks/lib/polymorphic_sites.rb region_3.aln >snps_region3.txt
      ruby ~/Code/Ruby/blocks/lib/polymorphic_sites.rb region_4.aln >snps_region4.txt
      ruby ~/Code/Ruby/blocks/lib/polymorphic_sites.rb region_5.aln >snps_region5.txt
      ruby ~/Code/Ruby/blocks/lib/polymorphic_sites.rb region_6.aln >snps_region6.txt

      echo "calculating snp densities for regions"
      bash ~/Code/Bash/snp_density.sh regions_positions.txt snps_region1.txt cluster_$clustr.aln >region1_snp_density.txt
      bash ~/Code/Bash/snp_density.sh regions_positions.txt snps_region2.txt cluster_$clustr.aln >region2_snp_density.txt
      bash ~/Code/Bash/snp_density.sh regions_positions.txt snps_region3.txt cluster_$clustr.aln >region3_snp_density.txt
      bash ~/Code/Bash/snp_density.sh regions_positions.txt snps_region4.txt cluster_$clustr.aln >region4_snp_density.txt
      bash ~/Code/Bash/snp_density.sh regions_positions.txt snps_region5.txt cluster_$clustr.aln >region5_snp_density.txt
      bash ~/Code/Bash/snp_density.sh regions_positions.txt snps_region6.txt cluster_$clustr.aln >region6_snp_density.txt

      #echo "listing mutation positions"
      #list the positions for each mutation
      #awk -F, '{print $1}' cluster_$clustr.mutations.txt | sort -nrk 1 |
      #uniq  >cluster_$clustr.mutation.positions.txt

      #echo "Sanitizing"
      #sed 's/A-G/G-A/g; s/T-C/C-T/g; s/A-T/T-A/g; s/A-C/C-A/g; s/C-G/G-C/g; s/T-G/G-T/g; s/--T/T--/g; s/--G/G--/g; s/--C/C--/g; s/--A/A--/g;' <cluster_$clustr.mutations.txt | 
      #awk -F, '{print $1,$2,$3,$4}' | sort -nrk 1 | uniq >cluster_$clustr.mutations.edited.txt

      #echo "Listing mutations per codon position"
      #awk '{print $2,$4}' <cluster_$clustr.mutations.edited.txt | sort | uniq -c >cluster_$clustr.mutations_percodon_position.txt

      echo "dN-dS substitutions"
       ~/Softwares/selection.rb -i cluster_$clustr.aln -S | awk '{print $4}' | sort | uniq -c >cluster_$clustr.dN_dS.txt
    cd ..
  done <most_clusters.txt

  echo "collect snps in $id identity threshold"
  #find . -type f -name *.mutations.edited.txt | xargs cat | awk '{print $2,$4}'| sort | uniq -c | sed -e 's/^[ \t]*//' >all.mutations.txt
  find . -type f -name *.snps.txt | xargs cat | sed -e 's/^[ \t]*//' >all.snps.txt

  echo "collect snps for regions"
  find . -type f -name snps_region1.txt | xargs cat | sed -e 's/^[ \t]*//' >all.region1.snps.txt
  find . -type f -name snps_region2.txt | xargs cat | sed -e 's/^[ \t]*//' >all.region2.snps.txt
  find . -type f -name snps_region3.txt | xargs cat | sed -e 's/^[ \t]*//' >all.region3.snps.txt
  find . -type f -name snps_region4.txt | xargs cat | sed -e 's/^[ \t]*//' >all.region4.snps.txt
  find . -type f -name snps_region5.txt | xargs cat | sed -e 's/^[ \t]*//' >all.region5.snps.txt

  echo "collect snp densities for regions"
  find . -type f -name region1_snp_density.txt | xargs cat | sed -e 's/^[ \t]*//' >all.region1.snps.density.txt
  find . -type f -name region2_snp_density.txt | xargs cat | sed -e 's/^[ \t]*//' >all.region2.snps.density.txt
  find . -type f -name region3_snp_density.txt | xargs cat | sed -e 's/^[ \t]*//' >all.region3.snps.density.txt
  find . -type f -name region4_snp_density.txt | xargs cat | sed -e 's/^[ \t]*//' >all.region4.snps.density.txt
  find . -type f -name region5_snp_density.txt | xargs cat | sed -e 's/^[ \t]*//' >all.region5.snps.density.txt
  find . -type f -name region6_snp_density.txt | xargs cat | sed -e 's/^[ \t]*//' >all.region6.snps.density.txt
  cd ..
done
