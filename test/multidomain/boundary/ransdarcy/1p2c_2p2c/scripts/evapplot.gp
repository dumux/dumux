reset
set term pngcairo size 800,600 solid
set output "./EvaporationRate.png"
set datafile separator ';'
DATA='./'
set ylabel "Evaporation Rate [mm/d]"
set xlabel "time [days]"
set yrange [0:7]
set xrange [0:8]
set key right top samplen 1

plot DATA.'storage_KEpsilon.csv'  u ($1/86400):4  w l lc rgb "red"    t 'KEpsilon',\
 DATA.'storage_LowReKEpsilon.csv' u ($1/86400):4  w l lc rgb "orange" t 'Low Re KEpsilon',\
 DATA.'storage_KOmega.csv'        u ($1/86400):4  w l lc rgb "green"  t 'KOmega',\
 DATA.'storage_OneEq.csv'         u ($1/86400):4  w l lc rgb "cyan"   t 'OneEq',\
 DATA.'storage_ZeroEq.csv'        u ($1/86400):4  w l lc rgb "purple" t 'ZeroEq'
