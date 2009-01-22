./richards ./grids/L_debug_pipe_small.dgf 30 10 > tmp1.txt
#grep "Sw-pC:" tmp1.txt | sed "s/Sw-pC://" > pC.csv
grep "pC-Sw:" tmp1.txt | sed "s/pC-Sw://" > Sw.csv
grep "pc-dSdpC:" tmp1.txt | sed "s/pc-dSdpC://" > dSdpC.csv 
#gnuplot 
#plot "pC.csv" with lines
#plot "Sw.csv" with lines
