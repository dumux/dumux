reset
set term pngcairo size 800,600 solid
set output "./VelocityDistribution_x4.png"
set datafile separator ','
DATA='./'
set xlabel "v_x/v_{x,max} [-]"
set ylabel "y [-]"
set yrange [0:5]
set xrange [-0.2:1.01]
set key left center samplen 1

plot  './references/bfs_experimentaldata_velocity.csv' u 9:8 w p t 'DriverSeegmiller Experimental Data' ,\
 './references/bfs_referencemodel_velocity.csv' u 9:8 w l lc rgb "black" t 'CFL3D Comparison Model'
