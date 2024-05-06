#!/bin/bash

#KEYS=`readlink -e $1`
#TSV=`readlink -e $2`
#BEST=`readlink -e $3`
#STATS=`readlink -e $4`
PYTHON=/home/ben/miniconda3/bin/python
KEYS=/home/ben/Projects/BCL_Workshop_2022/Tutorial_6/data/Runcie_chem_group_keys.txt
TSV=/home/ben/Projects/BCL_Workshop_2022/Tutorial_6/data/Runcie_et_al.data.mod.mapped.tsv
TSV_KD=/home/ben/Projects/BCL_Workshop_2022/Tutorial_6/data/Runcie_et_al.data.mod.mapped.kd.tsv
BEST=/home/ben/Projects/BCL_Workshop_2022/Tutorial_6/best.txt
STATS=/home/ben/Projects/BCL_Workshop_2022/Tutorial_6/scripts/basic_stats.py

# From base directory (directory-dependent)
for h in {L,V}; do rm -f best.L383${h}.txt; echo "key rosetta_ifx" >> best.L383${h}.txt ; for i in `cat $KEYS`; do grep ${i} $BEST | grep L383${h} ; done | awk '{print $1,$2}' >> best.L383${h}.txt ; done

rm -f V-L.dat ; echo "key V-L_rosetta_ifx" >> V-L.dat
paste best.L383L.txt best.L383V.txt | awk '{print $1,$4-$2}' | awk -F'/' '{print $1,$2}' | awk '{print $1,$3}' | tail -n+2 | sed "s:BRD2_L383L_::g" >> V-L.dat

# IC50
for i in `awk '{print $1}' $TSV` ; do rosetta=`grep ${i} V-L.dat | awk '{print $NF}'`; experiment=`grep ${i} $TSV | awk '{print $7}'`; echo $rosetta $experiment | awk '{if($2){print $0}}' ; done | grep -v "\-\.\-" | tail -n+2 | tr ' ' ',' > V-L.dat.csv

echo "IC50 stats: "
$PYTHON $STATS -input V-L.dat.csv -output data.stats -R -rho -mae -rmse -tau -plot -x_label "Experimental Binding Affinity (pKd)" -y_label "Predicted Interaction Energy (REU)"
cat data.stats
echo ""

# Kd
for i in `awk '{print $1}' $TSV_KD` ; do rosetta=`grep ${i} V-L.dat | awk '{print $NF}'`; experiment=`grep ${i} $TSV_KD | awk '{print $7}'`; echo $rosetta $experiment | awk '{if($2){print $0}}' ; done | grep -v "\-\.\-" | tail -n+2 | tr ' ' ',' > V-L.dat.kd.csv

echo "Kd stats: "
$PYTHON $STATS -input V-L.dat.kd.csv -output data.kd.stats -R -rho -mae -rmse -tau -plot -x_label "Experimental Binding Affinity (pKd)" -y_label "Predicted Interaction Energy (REU)"
cat data.kd.stats
echo ""
