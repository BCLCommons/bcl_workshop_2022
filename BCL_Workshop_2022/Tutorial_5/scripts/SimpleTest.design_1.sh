#!/bin/bash

####################################################################
#
# BCL-Rosetta sample RosettaScripts run file
#
# Author: Benjamin Brown
#
####################################################################

# Global variables
ROSETTA=/home/brownbp1/Rosetta/main/source/build/src/release/linux/3.10/64/x86/gcc/5.2/bcl/rosetta_scripts.bcl.linuxgccrelease
ROTAMER_LIBRARY=/hd1/brownbp1/workspace/bcl/rotamer_library
MUTABLE_FRAGMENTS=/home/brownbp1/BCL_Workshop_2022/Tutorial_5/inputs/ligands/LIG.sdf
MEDCHEM_FRAGMENTS=/home/brownbp1/BCL_Workshop_2022/Tutorial_5/inputs/medchem_fragments

# Input variables
XML=`readlink -e $1`
PROTEIN=`readlink -e $2`
LIGAND=`readlink -e $3`
PARAMS=`readlink -e $4`
PREFIX=$5

# Derived variables
protein=`basename $PROTEIN .pdb`
ligand=`basename $LIGAND .pdb`

# Run
$ROSETTA \
-parser:protocol $XML \
-in:file:s "$PROTEIN $LIGAND" \
-extra_res_fa ${PARAMS} \
-parser:script_vars rotamer_library="${ROTAMER_LIBRARY}" \
-parser:script_vars mutfrag="${MUTABLE_FRAGMENTS}" \
-parser:script_vars medchem_fragments="${MEDCHEM_FRAGMENTS}" \
-out:prefix $PREFIX \
-out:pdb_gz true \
-nstruct 100 \
-in:file:fullatom \
-restore_talaris_behavior \
-ignore_zero_occupancy false \
-linmem_ig 10 \
-constant_seed true \
-overwrite > ${PREFIX}.design_1.log
