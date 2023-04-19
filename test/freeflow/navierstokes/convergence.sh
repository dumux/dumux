#! /bin/bash
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# usage: build the executable in the build directory, copy the file to the folder with the executable, an example command to run the script is ./convergence.sh test_ff_stokes_donea ./params.input 2 4

if [ $# -ne 4 ]; then
  echo "usage: $0 EXECUTABLE_NAME PARAMS DIMENSION NUM_REFINEMENT_STEPS"
  exit
fi

make $1

LOGFILE=$1.log
rm $LOGFILE

DIM=$3
INITIALCELLS=4

if [ $DIM == 1 ]; then
GRIDCELLS="$INITIALCELLS"
elif [ $DIM == 2 ]; then
GRIDCELLS="$INITIALCELLS $INITIALCELLS"
else
GRIDCELLS="$INITIALCELLS $INITIALCELLS $INITIALCELLS"
fi

MAX=$4
MAXCELLS=$INITIALCELLS
for (( i=0; i <= $MAX; ++i )); do
  printf "refinement $i / $MAX "
  ((MAXCELLS *=2))
  ./$1 $2 -Grid.Cells "$GRIDCELLS" -Grid.Refinement $i -Problem.Name $1 -Problem.PrintConvergenceTestFile true
  echo "done."
done

echo "reset
set terminal pngcairo size 1000,750
set output 'convergence.png'
set xrange [($INITIALCELLS-1):$MAXCELLS]
set yrange [:1]
set format x '10^{%L}'
set format y '10^{%L}'
set xlabel 'number of cells per dimension'
set ylabel 'abs(l2 error)'
set log x
set log y" > $1.gp

PLOT="plot x**-1 w l lc rgb 'gray40' t 'order_1', \
x**-2 w l lc rgb 'blue20' t 'order_2', \
'$LOGFILE' u (\$4**(1./$DIM.)):10 w lp t 'pressure', \
'$LOGFILE' u (\$4**(1./$DIM.)):13 w lp t 'velocity_x'"
if [ $DIM == 2 ]; then
PLOT=$PLOT", '$LOGFILE' u (\$4**(1./$DIM.)):16 w lp t 'velocity_y'"
elif [ $DIM == 3 ]; then
PLOT=$PLOT", '$LOGFILE' u (\$4**(1./$DIM.)):16 w lp t 'velocity_y', '$LOGFILE' u (\$4**(1./$DIM.)):19 w lp t 'velocity_z'"
fi

echo $PLOT >> $1.gp
gnuplot --persist $1.gp
eog convergence.png
