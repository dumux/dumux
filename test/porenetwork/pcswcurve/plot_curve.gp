reset
set datafile missing "-nan"
set xlabel "S_w [-]"
set ylabel "p_c [P]"
set y2label "nThroats [-]"
set xrange [0:1.01]
set y2tics 50

plot './static_pc_sw_pc-s-curve.txt' w linespoints , './logfile_pcScurve.txt' u 3:6 w l, 'eqPoints_pcScurve.txt' u 3:6, './logfile_pcScurve.txt' u 3:2 w l, 'eqPoints_pcScurve.txt' u 3:7 axis x1y2, './static_pc_sw_pc-s-curve.txt' u 1:3 axis x1y2

pause mouse close
