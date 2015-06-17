reset
set datafile separator ';'

set xlabel 'Time [d]'
set ylabel 'Evaporation rate [mm/d]'
set xrange [0:5]
set yrange [0:5]
plot \
'storage.out' u ($1/86400):($3*86400) w l lw 2 t 'current'

set terminal pngcairo size 1200,900
set output 'evaporationRates.png'
replot