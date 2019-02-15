DATA="."
SCRIPT="../../../../../bin/postprocessing/extractlinedata.py"

runSim () {
./$1 $INPUT -Problem.Name $2 | tee -a logfile.out
input=`ls -ltr bfs*vtu | tail -n 1 | awk '{print $9}'`
echo $input" -> "$2 | tee -a logfile.out
pvpython $SCRIPT -f $input -o $DATA -of $2"_velocity_m4" -p1 $VProfile1_P1 -p2 $VProfile1_P2 -v 2 -r 10000 | tee -a logfile.out
pvpython $SCRIPT -f $input -o $DATA -of $2"_velocity_1" -p1 $VProfile1_P1 -p2 $VProfile1_P2 -v 2 -r 10000 | tee -a logfile.out
pvpython $SCRIPT -f $input -o $DATA -of $2"_velocity_4" -p1 $VProfile4_P1 -p2 $VProfile4_P2 -v 2 -r 10000 | tee -a logfile.out
pvpython $SCRIPT -f $input -o $DATA -of $2"_velocity_6" -p1 $VProfile6_P2 -p2 $VProfile6_P2 -v 2 -r 10000 | tee -a logfile.out
pvpython $SCRIPT -f $input -o $DATA -of $2"_velocity_10" -p1 $VProfile10_P2 -p2 $VProfile10_P2 -v 2 -r 10000 | tee -a logfile.out
pvpython $SCRIPT -f $input -o $DATA -of $2"_BaseProfile" -p1 $BaseProfile_P1 -p2 $BaseProfile_P2 -v 2 -r 10000 | tee -a logfile.out
}

### backwardsfacingstep
INPUT=scripts/params_verbose.input
VProfilem4_P1="-4.0 1.0 0.0"
VProfilem4_P2="-4.0 9.0 0.0"
VProfile1_P1="1.0 0.0 0.0"
VProfile1_P2="1.0 9.0 0.0"
VProfile4_P1="4.0 0.0 0.0"
VProfile4_P2="4.0 9.0 0.0"
VProfile6_P1="6.0 0.0 0.0"
VProfile6_P2="6.0 9.0 0.0"
VProfile10_P1="10.0 0.0 0.0"
VProfile10_P2="10.0 9.0 0.0"
BaseProfile_P1="0.0 0.0 0.0"
BaseProfile_P2="40.0 0.0 0.0"

runSim test_ff_rans_bfs_kepsilon      bfs_kepsilon
runSim test_ff_rans_bfs_komega        bfs_komega
runSim test_ff_rans_bfs_lowrekepsilon bfs_lowrekepsilon
runSim test_ff_rans_bfs_oneeq         bfs_oneeq
runSim test_ff_rans_bfs_zeroeq        bfs_zeroeq

gnuplot scripts/coefficientoffriction.gp
gnuplot scripts/velocitydistribution_m4.gp
gnuplot scripts/velocitydistribution_1.gp
gnuplot scripts/velocitydistribution_4.gp
gnuplot scripts/velocitydistribution_6.gp
gnuplot scripts/velocitydistribution_10.gp
eog *.png
