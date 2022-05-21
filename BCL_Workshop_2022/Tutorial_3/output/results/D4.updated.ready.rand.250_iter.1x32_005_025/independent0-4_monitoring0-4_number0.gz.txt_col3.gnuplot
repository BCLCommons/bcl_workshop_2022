set terminal pngcairo enhanced size 600,600
set output "./results/D4.updated.ready.rand.250_iter.1x32_005_025//independent0-4_monitoring0-4_number0.gz.txt_col3.data.FPR.png"
set origin 0.0,0.0
set nologscale
set grid
set autoscale fix
set key left bottom
set title "Contingency Matrix Measures"
set xlabel "FPR"
set ylabel ""
plot "./results/D4.updated.ready.rand.250_iter.1x32_005_025//independent0-4_monitoring0-4_number0.gz.txt_col3.data" using 7:1 with line title "InformationGainRatio" ls 1, \
  "./results/D4.updated.ready.rand.250_iter.1x32_005_025//independent0-4_monitoring0-4_number0.gz.txt_col3.data" using 7:2 with line title "TPR" ls 2, \
  "./results/D4.updated.ready.rand.250_iter.1x32_005_025//independent0-4_monitoring0-4_number0.gz.txt_col3.data" using 7:3 with line title "PPV" ls 3, \
  "./results/D4.updated.ready.rand.250_iter.1x32_005_025//independent0-4_monitoring0-4_number0.gz.txt_col3.data" using 7:4 with line title "LocalPPV" ls 4, \
  "./results/D4.updated.ready.rand.250_iter.1x32_005_025//independent0-4_monitoring0-4_number0.gz.txt_col3.data" using 7:6 with line title "Ideal-PPV" ls 6
reset
set terminal pngcairo enhanced size 600,600
set output "./results/D4.updated.ready.rand.250_iter.1x32_005_025//independent0-4_monitoring0-4_number0.gz.txt_col3.data.Cutoff.png"
set origin 0.0,0.0
set nologscale
set grid
set autoscale fix
set key left bottom
set title "Contingency Matrix Measures"
set xlabel "Cutoff"
set ylabel ""
plot "./results/D4.updated.ready.rand.250_iter.1x32_005_025//independent0-4_monitoring0-4_number0.gz.txt_col3.data" using 8:1 with line title "InformationGainRatio" ls 1, \
  "./results/D4.updated.ready.rand.250_iter.1x32_005_025//independent0-4_monitoring0-4_number0.gz.txt_col3.data" using 8:2 with line title "TPR" ls 2, \
  "./results/D4.updated.ready.rand.250_iter.1x32_005_025//independent0-4_monitoring0-4_number0.gz.txt_col3.data" using 8:3 with line title "PPV" ls 3, \
  "./results/D4.updated.ready.rand.250_iter.1x32_005_025//independent0-4_monitoring0-4_number0.gz.txt_col3.data" using 8:4 with line title "LocalPPV" ls 4, \
  "./results/D4.updated.ready.rand.250_iter.1x32_005_025//independent0-4_monitoring0-4_number0.gz.txt_col3.data" using 8:5 with line title "Ideal-PPV" ls 5
reset
set terminal pngcairo enhanced size 600,600
set output "./results/D4.updated.ready.rand.250_iter.1x32_005_025//independent0-4_monitoring0-4_number0.gz.txt_col3.data.logFPR.png"
set origin 0.0,0.0
set xrange [0.0001:1]
set logscale x
set grid
set autoscale fix
set key left bottom
set title "Contingency Matrix Measures"
set xlabel "FPR"
set ylabel ""
plot "./results/D4.updated.ready.rand.250_iter.1x32_005_025//independent0-4_monitoring0-4_number0.gz.txt_col3.data" using 7:1 with line title "InformationGainRatio" ls 1, \
  "./results/D4.updated.ready.rand.250_iter.1x32_005_025//independent0-4_monitoring0-4_number0.gz.txt_col3.data" using 7:2 with line title "TPR" ls 2, \
  "./results/D4.updated.ready.rand.250_iter.1x32_005_025//independent0-4_monitoring0-4_number0.gz.txt_col3.data" using 7:3 with line title "PPV" ls 3, \
  "./results/D4.updated.ready.rand.250_iter.1x32_005_025//independent0-4_monitoring0-4_number0.gz.txt_col3.data" using 7:4 with line title "LocalPPV" ls 4, \
  "./results/D4.updated.ready.rand.250_iter.1x32_005_025//independent0-4_monitoring0-4_number0.gz.txt_col3.data" using 7:6 with line title "Ideal-PPV" ls 6
reset
