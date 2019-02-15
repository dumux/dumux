reset
set term pngcairo size 800,600 solid
set output "./VelocityDistribution_x-4.png"
set datafile separator ','
DATA='./'
set xlabel "v_x/v_{x,max} [-]"
set ylabel "y [-]"
set yrange [1:5]
set xrange [0:1.025]
set key left center samplen 1

plot  './references/bfs_experimentaldata_velocity.csv' u 3:2 w p t 'DriverSeegmiller Experimental Data' ,\
 './references/bfs_referencemodel_velocity.csv' u 3:2 w l lc rgb "black" t 'CFL3D Comparison Model'
