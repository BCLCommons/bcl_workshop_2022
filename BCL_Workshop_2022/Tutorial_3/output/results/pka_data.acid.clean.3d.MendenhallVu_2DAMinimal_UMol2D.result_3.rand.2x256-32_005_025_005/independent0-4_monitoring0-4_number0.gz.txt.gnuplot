set terminal pngcairo enhanced size 500,500 # transparent size 2160,2160
set output "./results/pka_data.acid.clean.3d.MendenhallVu_2DAMinimal_UMol2D.result_3.rand.2x256-32_005_025_005//independent0-4_monitoring0-4_number0.gz.txt.gnuplot_txt_cross_correlation.png"
set encoding iso
set title "Cross Correlation"
set xlabel "Experimental"
set autoscale fix
plot x with lines linetype -1 linewidth 1 linecolor rgb 'black', \
  "./results/pka_data.acid.clean.3d.MendenhallVu_2DAMinimal_UMol2D.result_3.rand.2x256-32_005_025_005//independent0-4_monitoring0-4_number0.gz.txt.gnuplot_txt" using 1:2  with points pointsize 0.5 linecolor rgb 'red'
set title "Cross correlation: RMSD 1.69921 R^2 0.748258"