#!/bin/bash



# Global variables
ROSETTA=/home/ben/workspace/Rosetta/main/source/bin/rosetta_scripts.bcl.linuxgccrelease

# Input variables
XML=`readlink -e $1`
PROTEIN=`readlink -e $2`
LIGAND=`readlink -e $3`
PARAMS=`readlink -e $4`
RESFILE=`readlink -e $5`
PREFIX=$6
BUMP_R1=/home/ben/Projects/BCL_Workshop_2022/Tutorial_6/bump_medchem_fragments_r1/r1.all.sdf
BUMP_R2=/home/ben/Projects/BCL_Workshop_2022/Tutorial_6/bump_medchem_fragments_r2/r2.all.sdf

# Run
$ROSETTA \
-parser:protocol $XML \
-in:file:s "$PROTEIN $LIGAND" \
-parser:script_vars prefix="${PREFIX}" \
-parser:script_vars resfile=${RESFILE} \
-parser:script_vars methoxy_res=${METHOXY_RES} \
-parser:script_vars bump_r1=${BUMP_R1} \
-parser:script_vars bump_r2=${BUMP_R2} \
-parser:script_vars methoxy_scan_frag=/home/ben/Projects/BCL_Workshop_2022/Tutorial_6/methoxy_scan_fragment.sdf \
-parser:script_vars progress_file="${PREFIX}.gmc.log" \
-extra_res_fa ${PARAMS} \
-out:prefix $PREFIX \
-out:pdb_gz true \
-nstruct 10 \
-in:file:fullatom \
-restore_pre_talaris_2013_behavior \
-score:weights ligand \
-ignore_zero_occupancy false \
-linmem_ig 10 #> ${PREFIX}.log
