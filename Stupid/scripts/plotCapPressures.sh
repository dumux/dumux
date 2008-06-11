#! /bin/bash

BASEDIR=`dirname $0`

CSV_DIR=$BASEDIR/csv
SVG_DIR=$BASEDIR/svg
STUPID=$BASEDIR/../Stupid

if ! test -x "$STUPID"; then
    echo "Can't find an excutable 'Stupid' binary"
    exit 1;
fi

mkdir -p $CSV_DIR
mkdir -p $SVG_DIR

function vanGenuchen()
{
    ALPHA=$1
    N=$2
    OUT_FILE=$3
    

    $STUPID --useGenuchten \
        --genuchtenAlpha $ALPHA \
        --genuchtenN $N \
        --printCsvCapPressure \
   > $OUT_FILE
}

function brooksCorey()
{
    PE=$1
    ALPHA=$2
    OUT_FILE=$3
    

    $STUPID --useBrooksCorey \
        --brooksCoreyAlpha $ALPHA \
        --brooksCoreyPe $PE \
        --printCsvCapPressure \
   > $OUT_FILE
}


echo "Writing CSVs to $CSV_DIR"
N=3
GNUPLOT_GENUCHTEN_ALPHA=""
for ALPHA in 1 2 3 4 7 12; do
    OUT_FILE="$CSV_DIR/genuchten_a${ALPHA}_${N}.csv"
    vanGenuchen $ALPHA $N $OUT_FILE
    BLA="\"$OUT_FILE\" using 1:2 with lines title \"van Genuchten (alpha=$ALPHA, n=$N)\""
    if test -z "$GNUPLOT_GENUCHTEN_ALPHA"; then
        GNUPLOT_GENUCHTEN_ALPHA="$BLA"
    else
        GNUPLOT_GENUCHTEN_ALPHA="$GNUPLOT_GENUCHTEN_ALPHA, $BLA"
    fi
done

ALPHA=1.5
GNUPLOT_GENUCHTEN_N=""
for N in 1 2 3 4 7 12; do
    OUT_FILE="$CSV_DIR/genuchten_a${ALPHA}_${N}.csv"
    vanGenuchen $ALPHA $N $OUT_FILE
    BLA="\"$OUT_FILE\" using 1:2 with lines title \"van Genuchten (alpha=$ALPHA, n=$N)\""
    if test -z "$GNUPLOT_GENUCHTEN_N"; then
        GNUPLOT_GENUCHTEN_N="$BLA"
    else
        GNUPLOT_GENUCHTEN_N="$GNUPLOT_GENUCHTEN_N, $BLA"
    fi
done

PE=1
GNUPLOT_BROOKS_COREY_ALPHA=""
for ALPHA in 1 2 3 4 7 12; do
    OUT_FILE="$CSV_DIR/brookscorey_a${ALPHA}_Pe${N}.csv"
    brooksCorey $PE $ALPHA $OUT_FILE
    BLA="\"$OUT_FILE\" using 1:2 with lines title \"Books-Corey (alpha=$ALPHA, Pe=$PE)\""
    if test -z "$GNUPLOT_BROOKS_COREY_ALPHA"; then
        GNUPLOT_BROOKS_COREY_ALPHA="$BLA"
    else
        GNUPLOT_BROOKS_COREY_ALPHA="$GNUPLOT_BROOKS_COREY_ALPHA, $BLA"
    fi
done


echo "Plotting CSVs to the SVGs"
gnuplot <<EOF
set terminal svg
set output "$SVG_DIR/genuchten-alpha.svg"

# set logscale y
set yrange[0:4]

set xlabel "Saturation [%]"
set ylabel "p_c"

#set terminal png picsize 1024 768
#set output "resid-vs-iterations.pdf"

set title "Capillary Pressure vs. Effective Saturation "
plot $GNUPLOT_GENUCHTEN_ALPHA;
EOF

gnuplot <<EOF
set terminal svg
set output "$SVG_DIR/genuchten-n.svg"

# set logscale y
set yrange[0:4]

set xlabel "Saturation [%]"
set ylabel "p_c"

#set terminal png picsize 1024 768
#set output "resid-vs-iterations.pdf"

set title "Capillary Pressure vs. Effective Saturation "
plot $GNUPLOT_GENUCHTEN_N;
EOF

gnuplot <<EOF
set terminal svg
set output "$SVG_DIR/brooks-corey.svg"

# set logscale y
set yrange[0:4]

set xlabel "Saturation [%]"
set ylabel "p_c"

#set terminal png picsize 1024 768
#set output "resid-vs-iterations.pdf"

set title "Capillary Pressure vs. Effective Saturation "
plot $GNUPLOT_BROOKS_COREY_ALPHA;
EOF

