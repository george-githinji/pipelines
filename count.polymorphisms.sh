#!/bin/bash

function findPolymorps(){
  local dir=$1
  local change=$2
  
  find . -maxdepth 3 -type f -name "all.mutations.txt" |
  grep "/${dir}/" |
  xargs cat |
  grep -E "${change}" |
  awk '{sum+=$1} END {print sum}'
}

dir=(98 96 94 92 90 88)
for i in "${dir[@]}"
do
  ga1='G-A 1'
  ga2="G-A 2"
  ga3="G-A 3"

  gc1="G-C 1"
  gc2="G-C 2"
  gc3="G-C 3"

  gt1="G-T 1"
  gt2="G-T 2"
  gt3="G-T 3"

  ta1="T-A 1"
  ta2="T-A 2"
  ta3="T-A 3"

  ct1="C-T 1"
  ct2="C-T 2"
  ct3="C-T 3"

  ca1="C-A 1"
  ca2="C-A 2"
  ca3="C-A 3"

  echo "$(findPolymorps $i "$ga1")" "G-A" 1 $i
  echo "$(findPolymorps $i "$ga2")" "G-A" 2 $i
  echo "$(findPolymorps $i "$ga3")" "G-A" 3 $i
  
  echo "$(findPolymorps $i "$gc1")" "G-C" 1 $i
  echo "$(findPolymorps $i "$gc2")" "G-C" 2 $i
  echo "$(findPolymorps $i "$gc3")" "G-C" 3 $i
  
  echo "$(findPolymorps $i "$gt1")" "G-T" 1 $i
  echo "$(findPolymorps $i "$gt2")" "G-T" 2 $i
  echo "$(findPolymorps $i "$gt3")" "G-T" 3 $i
  
  echo "$(findPolymorps $i "$ca1")" "C-A" 1 $i
  echo "$(findPolymorps $i "$ca2")" "C-A" 2 $i
  echo "$(findPolymorps $i "$ca3")" "C-A" 3 $i
  
  echo "$(findPolymorps $i "$ta1")" "T-A" 1 $i
  echo "$(findPolymorps $i "$ta2")" "T-A" 2 $i
  echo "$(findPolymorps $i "$ta3")" "T-A" 3 $i
  
  echo "$(findPolymorps $i "$ct1")" "C-T" 1 $i
  echo "$(findPolymorps $i "$ct2")" "C-T" 2 $i
  echo "$(findPolymorps $i "$ct3")" "C-T" 3 $i
done
