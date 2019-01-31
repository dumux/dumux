reset
set term pngcairo size 800,600 solid
set output "./LawOfTheWall.png"
set datafile separator ','
DATA='./'
set xlabel "y^+ [-]"
set ylabel "u_+ [-]"
set yrange [0:25]
set log x
set xrange [1:1100]
set key right bottom samplen 1
set arrow from 5,0 to 5,25 lc rgb 'gray70' nohead
set arrow from 30,0 to 30,25 lc rgb 'gray70' nohead

plot './references/laufer_re50000_u+y+.csv' u 1:2 w p t 'Laufer 1954, Re=50000',\
 DATA.'lauferpipe_kepsilon.csv'      u 11:12 w l lc rgb "red"    t 'KEpsilon with u_{tau}',\
 DATA.'lauferpipe_kepsilon.csv'      u 15:16 w l lc rgb "orange" t 'KEpsilon with u_{tau,nom}',\
 DATA.'lauferpipe_lowrekepsilon.csv' u 11:12 w l lc rgb "yellow" t 'Low Re KEpsilon',\
 DATA.'lauferpipe_komega.csv'        u 11:12 w l lc rgb "green"  t 'KOmega',\
 DATA.'lauferpipe_oneeq.csv'         u 11:12 w l lc rgb "cyan"   t 'OneEq',\
 DATA.'lauferpipe_zeroeq.csv'        u 11:12 w l lc rgb "purple" t 'ZeroEq',\
  1/0.41*log(x)+5      w l lc rgb 'gray10',\
  x                    w l lc rgb 'gray40
