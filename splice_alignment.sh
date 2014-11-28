#!/bin/bash

positionsfile=$1
alignfile=$2

[ $# -eq 0 ] && { echo "Usage: $0 alignment file"; exit 1;}

counter=1
while read pos; do
    echo $pos
    cat $alignfile | extractalign -filter -regions $pos >region_$counter.aln
    counter=$(( $counter + 1))
done<$positionsfile


