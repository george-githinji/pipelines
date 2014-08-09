#!/bin/bash
db="/Users/george/454_data/blast-db/454"
classified="/Users/george/454_data/blast-db/454.classified.txt"

cluster=$1
position=$2
window=$3

clsname=$(basename "$cluster")
clsext="${clsname##*.}"
clsname="${clsname%.*}"

#remove suffix from alignment
~/Code/Ruby/Blocks/lib/extract regions --infile $cluster --pos $position --win $window | sort | \
  uniq | awk -v cl=$clsname.$position '{print ">"cl"." FNR "\n" $0}' \
  >"$clsname.$position.$window.fa"

blastn -query $clsname.$position.$window.fa -db $db -outfmt 6 |\
  awk -v win=$window '{if($3 == 100 && $4 == win) print $2}' | \
  sort -n | uniq >hits.454.$clsname.$position.$window.txt

grep -wFf hits.454.$clsname.$position.$window.txt $classified |\
  awk '{print $1,$7,$2}' | sort -k2 | sed 's/\..*$//g' | sort -n | uniq -c 
