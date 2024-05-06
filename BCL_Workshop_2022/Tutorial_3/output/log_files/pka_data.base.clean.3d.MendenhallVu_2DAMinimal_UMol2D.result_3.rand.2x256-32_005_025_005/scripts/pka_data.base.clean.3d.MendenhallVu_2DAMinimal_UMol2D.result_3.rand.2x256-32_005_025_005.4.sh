#!/bin/bash
cd /home/ben/Projects/pKa_prediction/first_pass
. ./log_files/pka_data.base.clean.3d.MendenhallVu_2DAMinimal_UMol2D.result_3.rand.2x256-32_005_025_005/scripts/functions.sh
PrologueStatus ./log_files/pka_data.base.clean.3d.MendenhallVu_2DAMinimal_UMol2D.result_3.rand.2x256-32_005_025_005/status.txt
run_command_check_status /home/ben/workspace/bcl/build/linux64_release/bin/bcl-apps-static.exe model:Train 'NeuralNetwork(balance=True,balance target ratio=0.1,balance max repeats=100000,transfer function = Rectifier(0.05),weight update = Simple(alpha=0.5,eta=0.0025),input dropout type=Zero,objective function =MAE_NMAD,input noise=0.0,iteration weight update=Attenuate(0.000,0,0,0.01),shuffle=True,steps per update=10,dropout(0.05,0.25,0.05),hidden architecture(256,32), rescale output dynamic range=True,rmsd report frequency=1,scaling=AveStd)'  -max_minutes 24000 -max_iterations 2000   -opencl Disable  --result_averaging_window 0 -random_seed `od -A n -t dI -N 4 /dev/urandom | tr -d -`  -final_objective_function 'MAE_NMAD'   -training 'Subset(number chunks=5,chunks="[0, 5) - [4] - [4]",filename="pka_data.base.clean.3d.MendenhallVu_2DAMinimal_UMol2D.result_3.rand.bin") '  -monitoring 'Subset(number chunks=5,chunks="[4]",filename="pka_data.base.clean.3d.MendenhallVu_2DAMinimal_UMol2D.result_3.rand.bin") '  -independent 'Subset(number chunks=5,chunks="[4]",filename="pka_data.base.clean.3d.MendenhallVu_2DAMinimal_UMol2D.result_3.rand.bin") '  --print_independent_predictions ./results/pka_data.base.clean.3d.MendenhallVu_2DAMinimal_UMol2D.result_3.rand.2x256-32_005_025_005/independent4_monitoring4_number0.gz  -storage_model 'File(directory=./models/pka_data.base.clean.3d.MendenhallVu_2DAMinimal_UMol2D.result_3.rand.2x256-32_005_025_005,prefix="model",write_descriptors=1,key=4)'  -logger File ./log_files/pka_data.base.clean.3d.MendenhallVu_2DAMinimal_UMol2D.result_3.rand.2x256-32_005_025_005/independent4_monitoring4_number0.txt &> /home/ben/Projects/pKa_prediction/first_pass/log_files/pka_data.base.clean.3d.MendenhallVu_2DAMinimal_UMol2D.result_3.rand.2x256-32_005_025_005/independent4_monitoring4_number0.job.txt