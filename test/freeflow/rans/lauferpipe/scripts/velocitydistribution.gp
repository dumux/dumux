reset
set term pngcairo size 800,600 solid
set output "./VelocityDistribution.png"
set datafile separator ','
DATA='./'
set xlabel "v_x/v_{x,max} [-]"
set ylabel "y [-]"
set yrange [0:0.5]
set key left center samplen 1

plot './references/laufer_re50000.csv' u 2:1 w p t 'Laufer 1954, Re=50000' ,\
 DATA.'lauferpipe_kepsilon.csv'      u 6:($29/0.2456) w l lc rgb "red"    t 'KEpsilon',\
 DATA.'lauferpipe_lowrekepsilon.csv' u 6:($25/0.2456) w l lc rgb "orange" t 'Low Re KEpsilon',\
 DATA.'lauferpipe_komega.csv'        u 6:($25/0.2456) w l lc rgb "green"  t 'KOmega',\
 DATA.'lauferpipe_oneeq.csv'         u 6:($24/0.2456) w l lc rgb "cyan"   t 'OneEq',\
 DATA.'lauferpipe_zeroeq.csv'        u 6:($23/0.2456) w l lc rgb "purple" t 'ZeroEq'
