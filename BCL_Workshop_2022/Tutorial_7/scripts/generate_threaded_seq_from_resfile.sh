#!/bin/bash

#################################
#
# Generates a sequence for threading with SimpleThreadingMover from a resfile
#
################################

N_RES=$1
RESFILE=`readlink -e $2`
OUTPUT=$3

for i in `seq 1 $N_RES`; do grep -w PIKAA $RESFILE | awk -v n="${i}" 'BEGIN{x="-"} {if(n==$1){x=$NF}} END {print x}' ; done | tr -d '\n' > $OUTPUT
