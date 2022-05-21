#!/bin/bash
#SBATCH --job-name=bcl-rosetta
#SBATCH --mem-per-cpu=7G
#SBATCH --time=12:00:00

# Global variables
ROSETTA=/home/ben/workspace/Rosetta/main/source/bin/rosetta_scripts.bcl.linuxgccrelease

# Input variables
XML=`readlink -e $1`
PROTEIN=`readlink -e $2`
LIGAND=`readlink -e $3`
PARAMS=`readlink -e $4`
RESFILE=`readlink -e $5`
BUMP_R1=`readlink -e $6`
BUMP_R2=`readlink -e $7`
METHOXY_RES=$8 # 16 or 17
HEAVY_PROB=$9
H_PROB=`python -c "print(1.0 - $HEAVY_PROB)"`
PREFIX=${10}
#PREFIX=${6}_${SLURM_ARRAY_TASK_ID}_

# Run
$ROSETTA \
-parser:protocol $XML \
-in:file:s "$PROTEIN $LIGAND" \
-parser:script_vars prefix="${PREFIX}" \
-parser:script_vars resfile=${RESFILE} \
-parser:script_vars methoxy_res=${METHOXY_RES} \
-parser:script_vars heavy_prob=${HEAVY_PROB} \
-parser:script_vars h_prob=${H_PROB} \
-parser:script_vars bump_r1=${BUMP_R1} \
-parser:script_vars bump_r2=${BUMP_R2} \
-parser:script_vars methoxy_scan_frag=/home/ben/Projects/BCL_Workshop_2022/Tutorial_6/methoxy_scan_fragment.sdf \
-extra_res_fa ${PARAMS} \
-out:prefix $PREFIX \
-out:pdb_gz true \
-nstruct 1 \
-in:file:fullatom \
-restore_pre_talaris_2013_behavior \
-score:weights ligand \
-ignore_zero_occupancy false \
-linmem_ig 10
