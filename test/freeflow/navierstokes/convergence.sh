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
set xlabel 'refinement'; \
set log y; \
set ylabel 'l2-error'; \
set arrow from 5,0.04 to 5,0.08 nohead lc 8; \
set arrow from 4,0.08 to 5,0.08 nohead lc 8; \
set arrow from 4,0.08 to 5,0.04 nohead lc 8" > $1.gp

USE_REL_ERROR=2 # set to "2" if the use relative error
PLOT="plot '$L2ERRORFILE' u :`expr 17 + $USE_REL_ERROR` w lp t 'pressure', '$L2ERRORFILE' u :`expr 23 + $USE_REL_ERROR` w lp t 'velocity'"
if [ $2 == 2 ]; then
PLOT=$PLOT", '$L2ERRORFILE' u :`expr 29 + $USE_REL_ERROR` w lp t 'velocity'"
elif [ $2 == 3 ]; then
PLOT=$PLOT", '$L2ERRORFILE' u :`expr 35 + $USE_REL_ERROR` w lp t 'velocity'"
fi

echo $PLOT >> $1.gp
gnuplot --persist $1.gp
