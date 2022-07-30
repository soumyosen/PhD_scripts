#!/bin/bash

for j in `seq 1 5`;
do

awk 'NR==38,NR==47 {print $2}' hyd$j.H.b.log >> hyd.50fr.H.b.log

done

