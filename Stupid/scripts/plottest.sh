#! /bin/bash

SNIPSNAP=`dirname $0`/snipsnap.py

rm *stream*_*.csv
cat test.csv | $SNIPSNAP

BLA=""
for TMP in *stream_*.csv; do
    if ! test -s "$TMP"; then
        continue
    fi

    if test -n "$BLA"; then
        BLA="$BLA, " 
    fi

    TITLE=`echo $TMP | sed "s/stream.*_\(.*\).csv/\1/"`
    BLA="$BLA \"$TMP\" using 1:2 with lines title \"$TITLE\""
done

gnuplot <<EOF
set terminal svg
set output  "test.svg"

#set logscale y
#set yrange[0:-0.003]
#set yrange[-1:1]
#set xrange[0.7:0.85]

set xlabel "Saturation [%]"
set ylabel "p_c"

#set terminal png picsize 1024 768
#set output "resid-vs-iterations.pdf"

set title "Capillary Pressure vs. Effective Saturation "
plot $BLA;
EOF

if test "$?" = "0"; then 
    inkscape test.svg &
fi
