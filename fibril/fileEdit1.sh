#!/bin/bash

#count=6

for j in `seq 1 5`;
do

awk 'NR==38,NR==47 {print $2}' hyd$j.H.b.log >> t.log

done

python data1.py 
#rm temp.log
rm t.log














