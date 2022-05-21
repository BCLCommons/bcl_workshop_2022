#!/bin/bash

#######################################
#
# Demo script to do structure-based design with BCL-Rosetta
# optimizing for target selectivity
#
# Author: Benjamin P. Brown
#
#######################################

# Global variables
ROSETTA=/home/ben/workspace/Rosetta/main/source/build/src/release/linux/5.10/64/x86/gcc/9/bcl/rosetta_scripts.bcl.linuxgccrelease
ROTAMER_LIBRARY=/home/ben/workspace/bcl/rotamer_library
MUTABLE_FRAGMENTS=/home/ben/Projects/BCL_Workshop_2022/Tutorial_5/inputs/ligands/LIG.sdf
MEDCHEM_FRAGMENTS=/home/ben/Projects/BCL_Workshop_2022/Tutorial_5/inputs/medchem_fragments
CONFS_LIST=/home/ben/Projects/BCL_Workshop_2022/Tutorial_5/inputs/receptor/conformers/confs.list
RINGS=/home/ben/Projects/BCL_Workshop_2022/Tutorial_7/input/rings/alternative_rings.west.sdf.gz

# Input variables
ROSETTA=/home/ben/workspace/Rosetta/main/source/bin/rosetta_scripts.bcl.linuxgccrelease
XML=`readlink -e $1`
PROTEIN=`readlink -e $2`
LIGAND=`readlink -e $3`
PARAMS=`readlink -e $4`
SEQ_NATIVE=`readlink -e $5`
SEQ2=`readlink -e $6`
SEQ3=`readlink -e $7`
SEQ5=`readlink -e $8`
PREFIX=$9

# Run
$ROSETTA \
-parser:protocol $XML \
-in:file:s "$PROTEIN $LIGAND" \
-parser:script_vars prefix="${PREFIX}" \
-parser:script_vars d4_seq=`cat ${SEQ_NATIVE}` \
-parser:script_vars d2_seq=`cat ${SEQ2}` \
-parser:script_vars d3_seq=`cat ${SEQ3}` \
-parser:script_vars d5_seq=`cat ${SEQ5}` \
-parser:script_vars rotamer_library="${ROTAMER_LIBRARY}" \
-parser:script_vars mutfrag="${MUTABLE_FRAGMENTS}" \
-parser:script_vars medchem_fragments="${MEDCHEM_FRAGMENTS}" \
-parser:script_vars rings="${RINGS}" \
-extra_res_fa ${PARAMS} \
-out:prefix $PREFIX \
-out:pdb false \
-out:pdb_gz false \
-nstruct 20 \
-in:file:fullatom \
-restore_talaris_behavior \
-ignore_zero_occupancy false \
-linmem_ig 10 \
-constant_seed \
-overwrite #> ${PREFIX}.log
