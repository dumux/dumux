# Test to see if the download and configuration worked properly

echo "DeprecationWarning: This script will be removed after 3.5!!!"
echo "The commands can be found in the handbook and are printed after installation"
echo "via installdumux.py"

cd dumux/dumux/build-cmake/test/porousmediumflow/1p/isothermal
make test_1p_tpfa
./test_1p_tpfa params.input
paraview *pvd
