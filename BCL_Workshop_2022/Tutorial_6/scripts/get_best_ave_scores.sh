#!/bin/bash

N=$1

# Go to base directory where all design directories are located
dir=$PWD
 
for i in `cat /home/ben/Projects/BCL_Workshop_2022/Tutorial_6/dir.list`; do
cd ${dir}/${i} ; 
ave_dg=`for j in \`ls *.pdb.gz\`; do dg=\`zgrep -w interface_delta_X $j | awk '{print $2}'\`; echo $i/$j $dg ; done | sort -gk2 | head -n ${N} | ${dir}/scripts/ave_column 2 | awk '{print $1}'` ;

colname=`for j in \`ls *.pdb.gz\`; do 
dg=\`zgrep -w interface_delta_X $j | awk '{print $2}'\`; 
echo $i/$j $dg ; 
done | sort -gk2 | head -n1 | awk '{print $1}'` ; 

echo $colname $ave_dg

done > ${dir}/best.txt
