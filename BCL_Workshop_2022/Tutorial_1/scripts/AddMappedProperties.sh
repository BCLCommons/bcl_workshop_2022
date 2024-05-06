#!/bin/sh
BCL=$1
INPUT=$2 # input sdf file
OUTPUT=$3 # output sdf file
MAP=$4 # must be in fixed-width .csv file; first column is key
NAME=$5 # property name
$BCL molecule:Properties -add 'Index' "Define(${NAME}=Mapped(delimiter=\",\", file=${MAP}, key=Cached(Index)))" "${NAME}" -input_filenames ${INPUT} -output ${OUTPUT}
