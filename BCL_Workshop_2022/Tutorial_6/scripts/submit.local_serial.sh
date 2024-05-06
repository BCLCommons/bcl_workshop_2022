#!/bin/bash

#############################################################
#							    #
# Wrapper submission script for BRD2 demo		    #
#							    #
# Run from /home/brownbp1/BRD/enumerative/		    #
#							    #
#############################################################

# Set current
dir=$PWD
ROUND=$1

if [[ $ROUND -eq 1 ]]; then
# 8-Methoxy Cl
for h in {V,L}; do for i in {0..3}; do for j in {0..2}; do k=1.0 ; l=16 ; mkdir ${dir}/BRD2_L383${h}_r1-${i}_r2-${j}_meo-8 ; cd ${dir}/BRD2_L383${h}_r1-${i}_r2-${j}_meo-8 ; bash ${dir}/scripts/ChemoGeneticsDesign.sh ${dir}/scripts/ChemoGeneticsDesign.xml ${dir}/BRD2.pdb ${dir}/XYZ.pdb ${dir}/XYZ.params ${dir}/resfiles/${h}.resfile ${dir}/bump_medchem_fragments_r1/r1_${i}.sdf ${dir}/bump_medchem_fragments_r2/r2_${j}.sdf ${l} ${k} ; done ; done ; done > /dev/null
fi

if [[ $ROUND -eq 2 ]]; then
# 9-Methoxy Cl
for h in {V,L}; do for i in {0..3}; do for j in {0..2}; do k=1.0 ; l=17 ; mkdir ${dir}/BRD2_L383${h}_r1-${i}_r2-${j}_meo-9 ; cd ${dir}/BRD2_L383${h}_r1-${i}_r2-${j}_meo-9 ; bash ${dir}/scripts/ChemoGeneticsDesign.sh ${dir}/scripts/ChemoGeneticsDesign.xml ${dir}/BRD2.pdb ${dir}/XYZ.pdb ${dir}/XYZ.params ${dir}/resfiles/${h}.resfile ${dir}/bump_medchem_fragments_r1/r1_${i}.sdf ${dir}/bump_medchem_fragments_r2/r2_${j}.sdf ${l} ${k} ; done ; done ; done > /dev/null
fi

if [[ $ROUND -eq 3 ]]; then
# 8-Methoxy Cl, bump1 is H
for h in {V,L}; do i=h ; for j in {0..2}; do k=0.0 ; l=16 ; mkdir ${dir}/BRD2_L383${h}_r1-${i}_r2-${j}_meo-8 ; cd ${dir}/BRD2_L383${h}_r1-${i}_r2-${j}_meo-8 ; bash ${dir}/scripts/ChemoGeneticsDesign.sh ${dir}/scripts/ChemoGeneticsDesign.xml ${dir}/BRD2.pdb ${dir}/XYZ.pdb ${dir}/XYZ.params ${dir}/resfiles/${h}.resfile ${dir}/bump_medchem_fragments_r1/r1_${i}.sdf ${dir}/bump_medchem_fragments_r2/r2_${j}.sdf ${l} ${k} ; done ; done > /dev/null
fi

if [[ $ROUND -eq 4 ]]; then
# 9-Methoxy Cl, bump1 is H
for h in {V,L}; do i=h ; for j in {0..2}; do k=0.0 ; l=17 ; mkdir ${dir}/BRD2_L383${h}_r1-${i}_r2-${j}_meo-9 ; cd ${dir}/BRD2_L383${h}_r1-${i}_r2-${j}_meo-9 ; bash ${dir}/scripts/ChemoGeneticsDesign.sh ${dir}/scripts/ChemoGeneticsDesign.xml ${dir}/BRD2.pdb ${dir}/XYZ.pdb ${dir}/XYZ.params ${dir}/resfiles/${h}.resfile ${dir}/bump_medchem_fragments_r1/r1_${i}.sdf ${dir}/bump_medchem_fragments_r2/r2_${j}.sdf ${l} ${k} ; done ; done > /dev/null
fi
