reset
set term pngcairo size 800,600 solid
set output "./FrictionCoefficient.png"
set datafile separator ','
DATA='./'
set xlabel "distance from the step [-]"
set ylabel "Coefficient of Friction at Base Wall"
set yrange [-0.15e-2:0.25e-2]
set xrange [0:40]
set key right bottom samplen 1

plot './references/bfs_experimentaldata_friction.csv' u 1:2 w p t 'DriverSeegmiller Experimental Data', \
 './references/bfs_referencemodel_friction.csv' u 1:2  w l lc rgb "black" t 'CFL3D Comparison Model',\
 DATA.'backwardfacingstep_friction_kepsilon.csv'      u 24:(($15+$16)/($5*0.5*44.2*44.2)) w l lc rgb "red"    t 'KEpsilon with u_{tau}',\
 DATA.'backwardfacingstep_friction_lowrekepsilon.csv' u 24:(($15+$16)/($5*0.5*44.2*44.2)) w l lc rgb "yellow" t 'Low Re KEpsilon',\
 DATA.'backwardfacingstep_friction_komega.csv'        u 24:(($15+$16)/($5*0.5*44.2*44.2)) w l lc rgb "green"  t 'KOmega',\
 DATA.'backwardfacingstep_friction_oneeq.csv'         u 24:(($15+$16)/($5*0.5*44.2*44.2)) w l lc rgb "cyan"   t 'OneEq',\
 DATA.'backwardfacingstep_friction_zeroeq.csv'        u 24:(($15+$16)/($5*0.5*44.2*44.2)) w l lc rgb "purple" t 'ZeroEq'
