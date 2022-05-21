#!/bin/bash
. ./log_files/D4.updated.ready.rand.250_iter.1x32_005_025/scripts/functions.sh
remove_duplicate_descriptor_files ./models/D4.updated.ready.rand.250_iter.1x32_005_025
run_merge_command /home/ben/workspace/bcl/build/linux64_release/bin/bcl-apps-static.exe model:PredictionMerge  -input_model_storage 'File(directory=./models/D4.updated.ready.rand.250_iter.1x32_005_025,prefix="model")'  -output ./results/D4.updated.ready.rand.250_iter.1x32_005_025/independent0-4_monitoring0-4_number0.gz.txt  -logger File ./log_files/D4.updated.ready.rand.250_iter.1x32_005_025/log_merge.txt  -opencl Disable 
run_merge_command /home/ben/workspace/bcl/build/linux64_release/bin/bcl-apps-static.exe model:ComputeStatistics  -message_level Verbose -input ./results/D4.updated.ready.rand.250_iter.1x32_005_025/independent0-4_monitoring0-4_number0.gz.txt -logger File ./log_files/D4.updated.ready.rand.250_iter.1x32_005_025/log_merge.txt  -opencl Disable  -filename_obj_function ./results/D4.updated.ready.rand.250_iter.1x32_005_025/final_objective.ind_merged.txt -obj_function 'AucRocCurve(cutoff=0.5,parity=1,x_axis_log=1,min fpr=0.001,max fpr=0.1)' 

