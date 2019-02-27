runSim () {
./$1 $INPUT -Vtk.OutputName $2 | tee -a logfile.out
}

### Evaporation Tests (Drying)
INPUT=scripts/params_evaporation.input

runSim test_md_boundary_darcy2p2cni_rans1p2cnikepsilon       evaporation_kepsilon
runSim test_md_boundary_darcy2p2cni_rans1p2cnikomega         evaporation_komega
runSim test_md_boundary_darcy2p2cni_rans1p2cnilowrekepsilon  evaporation_lowrekepsilon
runSim test_md_boundary_darcy2p2cni_rans1p2cnioneeq          evaporation_oneeq
runSim test_md_boundary_darcy2p2cni_rans1p2cnizeroeq         evaporation_zeroeq

gnuplot scripts/evapplot.gp
