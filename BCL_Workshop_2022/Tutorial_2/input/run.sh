# Features
for i in {dasatinib,amprenavir,ethinyl_estradiol}; do
	/home/ben/workspace/bcl/build/linux64_release/bin/bcl-apps-static.exe molecule:ConformerGenerator \
		-ensemble_filenames ${i}.sdf.gz \
		-conformers_single_file ../output/${i}.confs.sdf.gz \
		-max_iterations 8000 \
		-top_models 100 \
		-minimum_number_conformations 100 \
		-conformation_comparer RMSD 0.0 \
		-skip_cluster \
		-generate_3D
done

for i in {dasatinib,amprenavir,ethinyl_estradiol}; do 
	/home/ben/workspace/bcl/build/linux64_release/bin/bcl-apps-static.exe descriptor:GenerateDataset \
		-source "SdfFile(filename=../output/${i}.confs.sdf.gz)" \
		-feature_labels "Combine(3daSmoothSign(property=Atom_SigmaCharge,step size=0.25,temperature=100,steps=48,gaussian=1,interpolate=1))" \
		-result_labels "Constant(999)" \
		-output ../output/${i}.confs.3da.csv 
done

# Variance
for i in {dasatinib,amprenavir,ethinyl_estradiol}; do
	/home/ben/miniconda3/bin/python \
		/home/ben/workspace/pyben/var.py \
		-input ../output/${i}.confs.3da.csv \
		-output ../output/${i}.confs.3da.var.dat
done

# Plot
/home/ben/miniconda3/bin/python \
	../scripts/plot_variance.py \
	../output/dasatinib.confs.3da.var.dat \
	../output/amprenavir.confs.3da.var.dat \
	../output/ethinyl_estradiol.confs.3da.var.dat \
	../output/3da_variance.png

# Visualize
eog ../output/3da_variance.png

