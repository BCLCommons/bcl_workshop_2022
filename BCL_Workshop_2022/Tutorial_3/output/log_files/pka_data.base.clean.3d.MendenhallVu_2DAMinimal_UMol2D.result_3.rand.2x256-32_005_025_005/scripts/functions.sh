#!/bin/bash
export BCL_MAX_RERUNS=1
export BCL_CV_STATUS_FILENAME=./log_files/pka_data.base.clean.3d.MendenhallVu_2DAMinimal_UMol2D.result_3.rand.2x256-32_005_025_005/status.txt
export BCL_CV_STATUS_LOCK_FILENAME=/tmp/bcl_cv_9fc9b487e219aa186567071bb5785be9ef06f0a8.lock
export BCL_CV_WAKEUP_FILENAME=./log_files/pka_data.base.clean.3d.MendenhallVu_2DAMinimal_UMol2D.result_3.rand.2x256-32_005_025_005/scripts/wakeup
export BCL_CV_LOGFILES_PATH=./log_files/pka_data.base.clean.3d.MendenhallVu_2DAMinimal_UMol2D.result_3.rand.2x256-32_005_025_005
export PATH=/home/ben/Amber/amber20_src/bin::/usr/bin/:/usr/local/cuda-11.4/bin/:/home/ben/miniconda3/bin/:/usr/local/bin:/usr/sbin
export LD_LIBRARY_PATH=/home/ben/Amber/amber20_src/lib::/usr/local/cuda-11.4/lib64/:/usr/lib/:/usr/lib64/:/usr/lib32/:/usr/libx32/:/usr/libexec/
export PBS_O_WORKDIR=/home/ben/Projects/pKa_prediction/first_pass
export BCL_CV_SCRIPTS_PATH=./log_files/pka_data.base.clean.3d.MendenhallVu_2DAMinimal_UMol2D.result_3.rand.2x256-32_005_025_005/scripts
export FLOCK_MAX_WAIT=30
export BCL_MERGE_SCRIPT=./log_files/pka_data.base.clean.3d.MendenhallVu_2DAMinimal_UMol2D.result_3.rand.2x256-32_005_025_005/scripts/merge_ind_and_compute_results.sh
. /home/ben/workspace/bcl/scripts/machine_learning/tasks/util/functions.sh
