DATA="."
SCRIPT="../../../../bin/postprocessing/extractlinedata.py"

runSim () {
./$1 $INPUT -Problem.Name $2 | tee -a logfile.out
input=`ls -ltr lauferpipe*vtu | tail -n 1 | awk '{print $9}'`
echo $input" -> "$2 | tee -a logfile.out
pvpython $SCRIPT -f $input -o $DATA -of $2 -p1 $P1 -p2 $P2 -v 2 -r 10000 | tee -a logfile.out
}

### lauferpipe
INPUT=scripts/params_verbose.input
P1="8.0 0.0 0.0"
P2="8.0 0.12345 0.0"

runSim test_ff_rans_lauferpipe_kepsilon      lauferpipe_kepsilon
runSim test_ff_rans_lauferpipe_komega        lauferpipe_komega
runSim test_ff_rans_lauferpipe_lowrekepsilon lauferpipe_lowrekepsilon
runSim test_ff_rans_lauferpipe_oneeq         lauferpipe_oneeq
runSim test_ff_rans_lauferpipe_zeroeq        lauferpipe_zeroeq

gnuplot scripts/lawofthewall.gp
gnuplot scripts/velocitydistribution.gp
