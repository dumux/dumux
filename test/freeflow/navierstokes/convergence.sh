#! /bin/bash

if [ $# -ne 4 ]; then
  echo "usage: $0 EXECUTABLE_NAME PARAMS DIMENSION NUM_REFINEMENT_STEPS"
  exit
fi

make $1

LOGFILE=logfile_$1.out
rm $LOGFILE
L2ERRORFILE=l2error_$1.txt

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
  ./$1 params.input -Grid.Cells "$GRIDCELLS" -Grid.Refinement $i -Problem.PrintL2Error true &>> $LOGFILE
  echo "done."
done

grep "L2 error (abs/rel) for" $LOGFILE | tee $L2ERRORFILE
echo "reset
set terminal pngcairo size 1000,750
set output 'convergence.png'
set xrange [($INITIALCELLS-1):$MAXCELLS]
set yrange [:1]
set xlabel 'number of cells per dimension'
set ylabel 'l2 error (abs)'
set log x
set log y" > $1.gp

PLOT="plot x**-1 w l lc rgb 'gray40' t 'O_1', \
x**-2 w l lc rgb 'blue20' t 'O_2', \
'$L2ERRORFILE' u (\$6**(1./$DIM.)):17 w lp t 'pressure', \
'$L2ERRORFILE' u (\$6**(1./$DIM.)):23 w lp t 'velocity_x'"
if [ $DIM == 2 ]; then
PLOT=$PLOT", '$L2ERRORFILE' u (\$6**(1./$DIM.)):29 w lp t 'velocity_y'"
elif [ $DIM == 3 ]; then
PLOT=$PLOT", '$L2ERRORFILE' u (\$6**(1./$DIM.)):29 w lp t 'velocity_y', '$L2ERRORFILE' u (\$6**(1./$DIM.)):35 w lp t 'velocity_z'"
fi

echo $PLOT >> $1.gp
gnuplot --persist $1.gp
eog convergence.png
