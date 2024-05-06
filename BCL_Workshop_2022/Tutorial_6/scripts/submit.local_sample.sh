#!/bin/bash

#############################################################
#							    #
# Wrapper submission script for BRD2 benchmark		    #
#							    #
# Not robust; run from /home/brownbp1/BRD/enumerative/      #
#							    #
#############################################################

# Set current
dir=$PWD

# 8-Methoxy Cl
h=L; i=1; j=1; k=1.0 ; l=16 ; mkdir ${dir}/BRD2_L383${h}_r1-${i}_r2-${j}_meo-8 ; cd ${dir}/BRD2_L383${h}_r1-${i}_r2-${j}_meo-8 ; bash ${dir}/scripts/ChemoGeneticsDesign.sample.sh ${dir}/scripts/ChemoGeneticsDesign.xml ${dir}/BRD2.pdb ${dir}/XYZ.pdb ${dir}/XYZ.params ${dir}/resfiles/${h}.resfile ${dir}/bump_medchem_fragments_r1/r1_${i}.sdf ${dir}/bump_medchem_fragments_r2/r2_${j}.sdf ${l} ${k}
