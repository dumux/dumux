#! /bin/bash

if [ $# -ne 3 ]; then
  echo "usage: $0 EXECUTABLE_NAME DIMENSION NUM_REFINEMENT_STEPS"
  exit
fi

make $1

LOGFILE=logfile_$1.out
rm $LOGFILE
L2ERRORFILE=l2error_$1.txt

if [ $2 == 1 ]; then
GRIDCELLS="4"
elif [ $2 == 2 ]; then
GRIDCELLS="4 4"
else
GRIDCELLS="4 4 4"
fi

MAX=$3
for (( i=0; i <= $MAX; ++i )); do
  printf "refinement $i / $MAX "
  ./$1 $1.input -Grid.Cells "$GRIDCELLS" -Grid.Refinement $i -Problem.PrintL2Error true &>> $LOGFILE
  echo "done."
done

grep "L2 error (abs/rel) for" $LOGFILE | tee $L2ERRORFILE
echo "reset; \
set log x; \
set log y; \
set arrow from graph 0,1 to graph 1,0 nohead lc rgb 'gray'; \
set arrow from graph 0,1 to graph 1,0.5 nohead lc rgb 'gray'" > $1.gp

PLOT="plot '$L2ERRORFILE' u 6:17 w lp t 'pressure', '$L2ERRORFILE' u 6:23 w lp t 'velocity'"
if [ $2 == 2 ]; then
PLOT=$PLOT", '$L2ERRORFILE' u 6:29 w lp t 'velocity'"
elif [ $2 == 3 ]; then
PLOT=$PLOT", '$L2ERRORFILE' u 6:35 w lp t 'velocity'"
fi

echo $PLOT >> $1.gp
gnuplot --persist $1.gp
