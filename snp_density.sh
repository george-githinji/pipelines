#!/bin/bash
filename=$1  #regions.positions.txt
snpsfile=$2  #list of snps in the region
alignfile=$3 #alignment file

#extract the basename of the snps file
name=($basename $snpsfile)

#echo $name
#extract the region type from the basename of the file
n=$(echo $name | perl -pe '($_)=/(\d+)/')

#returns the size of the region under comparison
size=$(sed -ne "$n p" $filename | awk -F- '{print $2 - $1}')

#avoid indels | count the number of snps simply tally the them | remove the lead white space
totalsnps=$(awk '{if($4!="INDEL") print}' $snpsfile | awk 'END {print NR}' | sed -e 's/^[ \t]*//')

#get the number of sequences in the cluster alignment file
seqs=$(grep -c ">" $alignfile)
#echo "$seqs"

t=1
d=2
#calculate the number of pairwise comparisons N(N-1)/2
comp=$(($seqs*($seqs-$t)/$d))
#echo "$comparisons"

#calculate the density as the snps per region per number of sequences in the cluster
#density=$(bc<<< "scale=5;$totalsnps/$size*$seqs")  #this should take care of the depth of the snps

#pretty print with tab delimited columns
printf "%s\t%s\t%s\t%s\n" $size $totalsnps $seqs $comp
