#!/bin/bash
export BCL_MAX_RERUNS=1
export BCL_CV_STATUS_FILENAME=./log_files/D4.updated.ready.rand.250_iter.1x32_005_025/status.txt
export BCL_CV_STATUS_LOCK_FILENAME=/tmp/bcl_cv_d54bc47e7652e3d4922fcfd0b03b46d7858fccf0.lock
export BCL_CV_WAKEUP_FILENAME=./log_files/D4.updated.ready.rand.250_iter.1x32_005_025/scripts/wakeup
export BCL_CV_LOGFILES_PATH=./log_files/D4.updated.ready.rand.250_iter.1x32_005_025
export PATH=/home/ben/Amber/amber20_src/bin::/usr/bin/:/usr/local/cuda-11.4/bin/:/home/ben/miniconda3/bin/:/usr/local/bin:/usr/sbin
export LD_LIBRARY_PATH=/home/ben/Amber/amber20_src/lib::/usr/local/cuda-11.4/lib64/:/usr/lib/:/usr/lib64/:/usr/lib32/:/usr/libx32/:/usr/libexec/
export PBS_O_WORKDIR=/home/ben/Projects/BCL_Workshop_2022/Tutorial_3/output
export BCL_CV_SCRIPTS_PATH=./log_files/D4.updated.ready.rand.250_iter.1x32_005_025/scripts
export FLOCK_MAX_WAIT=30
export BCL_MERGE_SCRIPT=./log_files/D4.updated.ready.rand.250_iter.1x32_005_025/scripts/merge_ind_and_compute_results.sh
. /home/ben/workspace/bcl/scripts/machine_learning/tasks/util/functions.sh
