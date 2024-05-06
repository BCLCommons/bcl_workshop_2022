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
-out:prefix $PREFIX \
-out:pdb_gz true \
-nstruct 1 \
-in:file:fullatom \
-restore_talaris_behavior \
-ignore_zero_occupancy false \
-linmem_ig 10 \
-mute all \
-constant_seed true \
-overwrite
