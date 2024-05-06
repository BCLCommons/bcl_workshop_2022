#!/bin/bash
cd /net/tungsten/hd1/brownbp1/small_molecule_test_libraries/Mariusz_Benchmark_Set/confs
. ./log_files/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand/scripts/functions.sh
PrologueStatus ./log_files/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand/status.txt
run_command_check_status /net/tungsten/hd1/brownbp1/workspace/bcl/build/linux64_release/bin/bcl-apps-static.exe model:Train 'NeuralNetwork( transfer function = Sigmoid, weight update = Simple(alpha=0.5,eta=0.05),dropout(0.05,0.25),objective function = AucRocCurve(cutoff=0.5,parity=1,x_axis_log=1,min fpr=0.001,max fpr=0.1),scaling=AveStd,steps per update=1,hidden architecture(32),balance=True,balance target ratio=0.1,shuffle=True,input dropout type=Zero)'  -max_minutes 24000 -max_iterations 50   -opencl Disable  --result_averaging_window 0 -random_seed `od -A n -t dI -N 4 /dev/urandom | tr -d -`  -final_objective_function 'AucRocCurve(cutoff=0.5,parity=1,x_axis_log=1,min fpr=0.001,max fpr=0.1)'   -training 'Subset(number chunks=5,chunks="[0, 5) - [0] - [0]",filename="1834_combined.labeled.MendenhallMeiler2015.Minimal.rand.bin") '  -monitoring 'Subset(number chunks=5,chunks="[0]",filename="1834_combined.labeled.MendenhallMeiler2015.Minimal.rand.bin") '  -independent 'Subset(number chunks=5,chunks="[0]",filename="1834_combined.labeled.MendenhallMeiler2015.Minimal.rand.bin") '  --print_independent_predictions ./results/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand/independent0_monitoring0_number0.gz  -storage_model 'File(directory=./models/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand,prefix="model",write_descriptors=1,key=0)'  -logger File ./log_files/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand/independent0_monitoring0_number0.txt &> /net/tungsten/hd1/brownbp1/small_molecule_test_libraries/Mariusz_Benchmark_Set/confs/log_files/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand/independent0_monitoring0_number0.job.txt