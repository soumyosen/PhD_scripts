#!/bin/bash
count=4
#count2=10
#count3=10
#count4=10

while [ $count -lt 26 ]
do
        cd dist$count

        cat r2.colvars.traj  > us.colvars.traj

        rm -rf dist.txt
        awk '(NR>1) && (($0 !~ /step/)) {print $2}' us.colvars.traj > dist.txt
#        awk '{print $5}' sample_rna.job0.$count.sort.history > rc.sort
        cp ../hist.py .
        python3 hist.py
        mv dist.hist dist$count.hist
        cp dist$count.hist ..
        cd ..

     count=`expr $count + 1`
done


#while [ $count2 -lt 12 ]
#do
#        cd zz-$count2
#
#        cat us-0.colvars.traj  > us.colvars.traj
#
#        rm -rf dist.txt
#        awk '(NR>1) && (($0 !~ /step/)) {print $2}' us.colvars.traj > dist.txt
##        awk '{print $5}' sample_rna.job0.$count.sort.history > rc.sort
#        cp ../hist.py .
#        python3 hist.py
#
#        cd ..
#
#     count2=`expr $count2 + 1`
#done

#while [ $count3 -lt 16 ]
#do
#        cd zz-$count3
#
#	cat z-0.colvars.traj z-1.colvars.traj > z.colvars.traj
#
#        rm -rf dist.txt
#        awk '(NR>1) && (($0 !~ /step/)) {print $2}' z.colvars.traj > dist.txt
#        #awk '{print $5}' sample_rna.job0.$count.sort.history > rc.sort
#        cp ../hist.py .
#        python3 hist.py
#
#        cd ..
#
#     count3=`expr $count3 + 1`
#done
#
#
#while [ $count4 -lt 14 ]
#do
#        cd zz-n$count4
#
#	cat z-0.colvars.traj z-1.colvars.traj > z.colvars.traj
#
#        rm -rf dist.txt
#        awk '(NR>1) && (($0 !~ /step/)) {print $2}' z.colvars.traj > dist.txt
#        #awk '{print $5}' sample_rna.job0.$count.sort.history > rc.sort
#        cp ../hist.py .
#        python3 hist.py
#
#        cd ..
#
#     count4=`expr $count4 + 1`
#done
#
