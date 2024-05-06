#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=10:00:00

export LD_LIBRARY_PATH="/dors/meilerlab/apps/Linux2/x86_64/lib64:/dors/meilerlab/apps/Linux2/x86_64/lib:/net/antimony/hd0/vuot2/Rosetta/Rosetta/main/source/build/external/release/linux/3.10/64/x86/gcc/5.2/default/"

#parameters
PREFIX=$1
TEMPLATE=`readlink -e $2`
XML=`readlink -e $3`
RESFILE=`readlink -e $4`
RESS=$5 #list of ncaa names seperated by ","
PARAMS=$6 #list of ncaa params files seperated by " "
NSTRUCTS=$7
OUTDIR=`readlink -e $8`

#Rosetta version
ROSETTA=/mnt/ssd1/Rosetta/main/source/bin/rosetta_scripts.bcl.linuxgccrelease
#ROSETTA=/ors/meilerlab/apps/rosetta/rosetta-3.12/main/source/bin/rosetta_scripts.linuxgccrelease
#ROSETTA=/path/to/rosetta/main/source/bin/rosetta_scripts.linuxgccrelease
#run rosetta
cd $OUTDIR
mkdir -p output_files_${PREFIX}/
# Run design with the template 
$ROSETTA -parser:protocol $XML \
    -out:pdb_gz \
    -out:prefix ${PREFIX} \
    -nstruct $NSTRUCTS \
    -out:file:scorefile ${PREFIX}_scores.out \
    -s $TEMPLATE \
    -native $TEMPLATE \
    -out:path:all ${OUTDIR}/output_files_${PREFIX}/ \
    -ex1 \
    -ex2 \
    -extra_res_fa ${PARAMS}\
    -parser:script_vars res_type=${RESS} resf_file=${RESFILE}
