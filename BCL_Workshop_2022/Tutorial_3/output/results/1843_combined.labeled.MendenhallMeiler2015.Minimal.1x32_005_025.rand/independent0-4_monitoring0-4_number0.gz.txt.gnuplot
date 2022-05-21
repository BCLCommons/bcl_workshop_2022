set terminal pngcairo enhanced size 600,600
set output "./results/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand//independent0-4_monitoring0-4_number0.gz.txt.data.FPR.png"
set origin 0.0,0.0
set nologscale
set grid
set autoscale fix
set key left bottom
set title "Contingency Matrix Measures"
set xlabel "FPR"
set ylabel ""
plot "./results/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand//independent0-4_monitoring0-4_number0.gz.txt.data" using 7:1 with line title "InformationGainRatio" ls 1, \
  "./results/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand//independent0-4_monitoring0-4_number0.gz.txt.data" using 7:2 with line title "TPR" ls 2, \
  "./results/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand//independent0-4_monitoring0-4_number0.gz.txt.data" using 7:3 with line title "PPV" ls 3, \
  "./results/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand//independent0-4_monitoring0-4_number0.gz.txt.data" using 7:4 with line title "LocalPPV" ls 4, \
  "./results/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand//independent0-4_monitoring0-4_number0.gz.txt.data" using 7:6 with line title "Ideal-PPV" ls 6
reset
set terminal pngcairo enhanced size 600,600
set output "./results/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand//independent0-4_monitoring0-4_number0.gz.txt.data.Cutoff.png"
set origin 0.0,0.0
set nologscale
set grid
set autoscale fix
set key left bottom
set title "Contingency Matrix Measures"
set xlabel "Cutoff"
set ylabel ""
plot "./results/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand//independent0-4_monitoring0-4_number0.gz.txt.data" using 8:1 with line title "InformationGainRatio" ls 1, \
  "./results/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand//independent0-4_monitoring0-4_number0.gz.txt.data" using 8:2 with line title "TPR" ls 2, \
  "./results/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand//independent0-4_monitoring0-4_number0.gz.txt.data" using 8:3 with line title "PPV" ls 3, \
  "./results/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand//independent0-4_monitoring0-4_number0.gz.txt.data" using 8:4 with line title "LocalPPV" ls 4, \
  "./results/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand//independent0-4_monitoring0-4_number0.gz.txt.data" using 8:5 with line title "Ideal-PPV" ls 5
reset
set terminal pngcairo enhanced size 600,600
set output "./results/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand//independent0-4_monitoring0-4_number0.gz.txt.data.logFPR.png"
set origin 0.0,0.0
set xrange [0.0001:1]
set logscale x
set grid
set autoscale fix
set key left bottom
set title "Contingency Matrix Measures"
set xlabel "FPR"
set ylabel ""
plot "./results/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand//independent0-4_monitoring0-4_number0.gz.txt.data" using 7:1 with line title "InformationGainRatio" ls 1, \
  "./results/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand//independent0-4_monitoring0-4_number0.gz.txt.data" using 7:2 with line title "TPR" ls 2, \
  "./results/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand//independent0-4_monitoring0-4_number0.gz.txt.data" using 7:3 with line title "PPV" ls 3, \
  "./results/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand//independent0-4_monitoring0-4_number0.gz.txt.data" using 7:4 with line title "LocalPPV" ls 4, \
  "./results/1834_combined.labeled.MendenhallMeiler2015.Minimal.1x32_005_025.rand//independent0-4_monitoring0-4_number0.gz.txt.data" using 7:6 with line title "Ideal-PPV" ls 6
reset
