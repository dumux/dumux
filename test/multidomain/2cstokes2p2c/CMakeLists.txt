# Test for 2cstokes2p2c can be performed by the following commands:
# make test_2cstokes2p2c && ./test_2cstokes2p2c -ParameterFile test_references/test_2cstokes2p2c_reference.input | tee log_test.out
# ../../../bin/runTest.sh ../../../bin/fuzzycomparevtu.py test_references/stokes2c-00043.vtu stokes2c-00043.vtu /usr/bin/echo " "
# ../../../bin/runTest.sh ../../../bin/fuzzycomparevtu.py test_references/2p2c-00043.vtu 2p2c-00043.vtu /usr/bin/echo " "
add_dumux_test(test_2cstokes2p2c test_2cstokes2p2c test_2cstokes2p2c
  ${CMAKE_CURRENT_BINARY_DIR}/test_2cstokes2p2c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_2cstokes2p2c.input)
