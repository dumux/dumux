#! /bin/bash

MYDIR=`dirname $0`
SNIPSNAP=$MYDIR/snipsnap.py
EXPDATADIR=$MYDIR/../Problems/Lenhard/exp1
SIMNAME=LenhardExp1
SIMOUTPUT=$SIMNAME.log

function plotAtPos()
{
    OUTFILENAME=$1
    TITLE=$2
    FILE_EXP=$3
    FILE_SIM=$4

    BLA=""
    BLA="\"$FILE_EXP\" using 1:2 with dots title \"Experiment\""
    BLA="$BLA, \"$FILE_SIM\" using 1:2 with lines title \"Simulation\""

    gnuplot <<EOF
set terminal svg
set output  "$OUTFILENAME"

set xrange[0:75.0]
set yrange[0:1]

set title "$TITLE"
set xlabel "Time [h]"
set ylabel "Water saturation [%]"

plot $BLA;
EOF
}

# run simulation 
echo "Running simulation.."
cmake -DCMAKE_BUILD_TYPE=release .
make $SIMNAME
./Problems/Lenhard/$SIMNAME 1e100 10 | tee $SIMOUTPUT | grep -Ei "timestep"

echo "Separating output.."
rm -f stream*_*.csv
cat $SIMOUTPUT | $SNIPSNAP > /dev/null

echo "Generating SVG files.."
plotAtPos "exp1_67cm.svg" "Experiment I, pos=67cm" "$EXPDATADIR/67cm_exp1.csv" "stream_pos=0.67.csv"
plotAtPos "exp1_57cm.svg" "Experiment I, pos=57cm" "$EXPDATADIR/57cm_exp1.csv" "stream_pos=0.57.csv"
plotAtPos "exp1_47cm.svg" "Experiment I, pos=47cm" "$EXPDATADIR/47cm_exp1.csv" "stream_pos=0.47.csv"
plotAtPos "exp1_37cm.svg" "Experiment I, pos=37cm" "$EXPDATADIR/37cm_exp1.csv" "stream_pos=0.37.csv"
plotAtPos "exp1_27cm.svg" "Experiment I, pos=27cm" "$EXPDATADIR/27cm_exp1.csv" "stream_pos=0.27.csv"

rm *.csv
