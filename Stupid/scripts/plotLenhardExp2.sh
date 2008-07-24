#! /bin/bash

MYDIR=`dirname $0`
SNIPSNAP=$MYDIR/snipsnap.py
EXPDATADIR=$MYDIR/../Problems/Lenhard/exp2
SIMNAME=LenhardExp2
SIMOUTPUT=$SIMNAME.log

function plotAtPos()
{
    OUTFILENAME=$1
    TITLE=$2
    FILE_EXP=$3
    FILE_SIM=$4

    SWAPPSIM_CSV=`dirname $FILE_SIM`/VgSwapp_`basename $FILE_SIM`
    BLA=""
    BLA="\"$FILE_SIM\" using 1:2 with lines title \"PL Simulation \""
    if test -r $SWAPPSIM_CSV; then
        BLA="$BLA, \"$SWAPPSIM_CSV\" using 1:2 with lines title \"VG Simulation w. Snr#\""
    fi

    if test -n "$FILE_EXP"; then
        BLA="$BLA, \"$FILE_EXP\" using 1:2 with dots title \"Experiment\""
    fi
    
    gnuplot <<EOF
set terminal svg
set output  "$OUTFILENAME"

$XRANGE
$YRANGE

set title "$TITLE"
set xlabel "$XLABEL"
set ylabel "$YLABEL"

plot $BLA;
EOF
}

# run simulation 
echo "Running simulation.."
cmake -DCMAKE_BUILD_TYPE=release .
make $SIMNAME
./Problems/Lenhard/$SIMNAME 1e100 10 | tee $SIMOUTPUT | grep -Ei "timestep"

echo "Separating output.."
#rm -f stream*_*.csv
rm exp2_*.svg
#rm *stream_*.csv
cat $SIMOUTPUT | $SNIPSNAP > /dev/null

echo "Generating SVG files.."
XRANGE="set xrange[0:10.0]"
YRANGE="set yrange[0:1]"
XLABEL="Time [h]"
YLABEL="Water saturation"
plotAtPos "exp2_70cm.svg" "Experiment II, pos=70cm" "$EXPDATADIR/70cm_exp2.csv" *"stream_pos=0.70.csv"
plotAtPos "exp2_60cm.svg" "Experiment II, pos=60cm" "$EXPDATADIR/60cm_exp2.csv" *"stream_pos=0.60.csv"
plotAtPos "exp2_50cm.svg" "Experiment II, pos=50cm" "$EXPDATADIR/50cm_exp2.csv" *"stream_pos=0.50.csv"
plotAtPos "exp2_40cm.svg" "Experiment II, pos=40cm" "$EXPDATADIR/40cm_exp2.csv" *"stream_pos=0.40.csv"
plotAtPos "exp2_30cm.svg" "Experiment II, pos=30cm" "$EXPDATADIR/30cm_exp2.csv" *"stream_pos=0.30.csv"

XRANGE="set xrange[0:1.0]"
YRANGE=""
XLABEL="Water saturation"
YLABEL="Capillary pressure [Pa]"
plotAtPos "exp2_pC_70cm.svg" "Experiment II, pC at 70cm" "$EXPDATADIR/pC_70cm_exp2.csv" *"stream_pC_pos=0.70.csv" 
plotAtPos "exp2_pC_60cm.svg" "Experiment II, pC at 60cm" "$EXPDATADIR/pC_60cm_exp2.csv" *"stream_pC_pos=0.60.csv" 
plotAtPos "exp2_pC_50cm.svg" "Experiment II, pC at 50cm" "$EXPDATADIR/pC_50cm_exp2.csv" *"stream_pC_pos=0.50.csv" 
plotAtPos "exp2_pC_40cm.svg" "Experiment II, pC at 40cm" "$EXPDATADIR/pC_40cm_exp2.csv" *"stream_pC_pos=0.40.csv"

