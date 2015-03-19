#!/bin/bash
while read line 
do
  #echo $line 
  awk -v r=$line '{split($0,a,"."); print a[1],$1,r}' < cluster_$line/members.txt
done < $1

