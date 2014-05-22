# Test for 2cnistokes2p2cni can be performed by the following commands:
# make test_2cnistokes2p2cni && ./test_2cnistokes2p2cni -ParameterFile test_references/test_2cnistokes2p2cni_reference.input | tee log_test.out
# ../../../bin/runTest.sh ../../../bin/fuzzycomparevtu.py test_references/stokes2cni-00043.vtu stokes2cni-00043.vtu /usr/bin/echo " "
# ../../../bin/runTest.sh ../../../bin/fuzzycomparevtu.py test_references/2p2cni-00043.vtu 2p2cni-00043.vtu /usr/bin/echo " "
add_dumux_test(test_2cnistokes2p2cni test_2cnistokes2p2cni test_2cnistokes2p2cni
  ${CMAKE_CURRENT_BINARY_DIR}/test_2cnistokes2p2cni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_2cnistokes2p2cni.input)
