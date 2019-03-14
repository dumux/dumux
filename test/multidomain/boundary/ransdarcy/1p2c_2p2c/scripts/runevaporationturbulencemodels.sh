runSim () {
./$1 $INPUT -Vtk.OutputName $2 -RANS.Problem.InletVelocity $3| tee -a logfile.out
}

### Evaporation Tests (Drying)
INPUT=scripts/params_evaporation.input

runSim test_md_boundary_darcy2p2cni_rans1p2cnikepsilon       evaporation_kepsilon_3       3.0
runSim test_md_boundary_darcy2p2cni_rans1p2cnikomega         evaporation_komega_3         3.0
runSim test_md_boundary_darcy2p2cni_rans1p2cnilowrekepsilon  evaporation_lowrekepsilon_3  3.0
runSim test_md_boundary_darcy2p2cni_rans1p2cnioneeq          evaporation_oneeq_3          3.0
runSim test_md_boundary_darcy2p2cni_rans1p2cnizeroeq         evaporation_zeroeq_3         3.0
runSim test_md_boundary_darcy2p2cni_rans1p2cnikomega         evaporation_komega_0p3       0.3
runSim test_md_boundary_darcy2p2cni_rans1p2cnikomega         evaporation_komega_15        15.0
runSim test_md_boundary_darcy2p2cni_rans1p2cnikomega         evaporation_komega_30        30.0


gnuplot scripts/plotevaporationrates.gp
