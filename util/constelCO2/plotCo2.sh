#!/bin/bash

./printcsv > tmp.csv

gnuplot & <<EOF
#set terminal svg
#set output "$SVG_DIR/genuchten-n.svg"

# set logscale y
#set yrange[0:4]

set xlabel "T [K]"
set ylabel "p [Pa]"
set ylabel "z [rho|sol|enth]"

#set terminal png picsize 1024 768
#set output "resid-vs-iterations.pdf"

set title "Constrel CO2 "
set style data lines
splot "tmp.csv"
EOF


