#!/bin/bash

####################################################################
#
# BCL-Rosetta sample RosettaScripts run file
#
# Author: Benjamin Brown
#
####################################################################

# Global variables
ROSETTA=/home/ben/workspace/Rosetta/main/source/build/src/release/linux/5.10/64/x86/gcc/9/bcl/rosetta_scripts.bcl.linuxgccrelease
ROTAMER_LIBRARY=/home/ben/workspace/bcl/rotamer_library
MUTABLE_FRAGMENTS=/home/ben/Projects/BCL_Workshop_2022/Tutorial_5/inputs/ligands/LIG.sdf
MEDCHEM_FRAGMENTS=/home/ben/Projects/BCL_Workshop_2022/Tutorial_5/inputs/medchem_fragments
CONFS_LIST=/home/ben/Projects/BCL_Workshop_2022/Tutorial_5/inputs/receptor/conformers/confs.list

# Input variables
XML=`readlink -e $1`
PROTEIN=`readlink -e $2`
LIGAND=`readlink -e $3`
PARAMS=$4
TEMP=$5
TYPE=$6
PREFIX=$7

# Derived variables
protein=`basename $PROTEIN .pdb`
ligand=`basename $LIGAND .pdb`

# Run
$ROSETTA \
-parser:protocol $XML \
-in:file:s "$PROTEIN $LIGAND" \
-extra_res_fa "$PARAMS".fa.params \
-extra_res_cen "$PARAMS".cen.params \
-parser:script_vars rotamer_library="${ROTAMER_LIBRARY}" \
-parser:script_vars mutfrag="${MUTABLE_FRAGMENTS}" \
-parser:script_vars medchem_fragments="${MEDCHEM_FRAGMENTS}" \
-parser:script_vars temp="${TEMP}" \
-parser:script_vars tki_type="${TYPE}" \
-parser:script_vars confs_list="${CONFS_LIST}" \
-out:prefix $PREFIX \
-out:pdb_gz true \
-nstruct 10 \
-in:file:fullatom \
-restore_talaris_behavior \
-ignore_zero_occupancy false \
-linmem_ig 10 \
-mute protocols.qsar.scoring_grid.GridManager \
-constant_seed true \
-overwrite > ${PREFIX}.log
