#!/bin/bash

count=6

for j in `seq 1 5`;
do
awk '/^B/' hyd$j.BH.c1.log > temp.log

for i in `seq 0 9`;
do
let a=" ($count*$i)+1 "
let b=" $count*($i+1) "
let c=" $i+1 "
#echo "$a $b"
awk -v "linestart=$a" -v "lineend=$b" 'NR==linestart,NR==lineend {print $3}' temp.log >> t.log
done

done

python data.py 
rm temp.log
rm t.log














