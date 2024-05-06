#!/bin/bash

# Go to base directory where all design directories are located
dir=$PWD

for i in `ls -d BRD2_*`; do cd ${dir}/${i} ; for j in `ls *.pdb.gz`; do dg=`zgrep -w interface_delta_X $j | awk '{print $2}'`; echo $i/$j $dg ; done | sort -gk2 | head -n1 ; done > ${dir}/best.txt
#for i in `cat /home/ben/BRD/dir.list`; do cd ${dir}/${i} ; for j in `ls *.pdb.gz`; do dg=`zgrep -w interface_delta_X $j | awk '{print $2}'`; echo $i/$j $dg ; done | sort -gk2 | head -n1 ; done > ${dir}/best.txt
