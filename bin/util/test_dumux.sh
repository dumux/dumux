# Test to see if the download and configuration worked properly

cd DUMUX/dumux/build-cmake/test/porousmediumflow/1p/implicit/isothermal
make test_1p_tpfa
./test_1p_tpfa params.input
paraview *pvd
