#!/bin/bash

# Global variables
ROSETTA=/home/ben/workspace/Rosetta/main/source/build/src/release/linux/5.10/64/x86/gcc/9/bcl/rosetta_scripts.bcl.linuxgccrelease

# Input variables
XML=`readlink -e $1`
INPUT=`readlink -e $2`
PARAMS=`readlink -e $3`
POS=$4
RES=$5
PREFIX=$6

# Run
$ROSETTA \
-parser:protocol $XML \
-in:file:s "$INPUT" \
-parser:script_vars prefix="${PREFIX}" \
-parser:script_vars pos="${POS}" \
-parser:script_vars res="${RES}" \
-extra_res_fa ${PARAMS} \
-out:prefix $PREFIX \
-use_input_sc true \
-out:pdb_gz true \
-nstruct 1 \
-in:file:fullatom \
-ignore_zero_occupancy false \
-linmem_ig 10 \
-in:detect_disulf \
-relax:constrain_relax_to_start_coords \
-mute core.select.residue_selector.PrimarySequenceNeighborhoodSelector \
