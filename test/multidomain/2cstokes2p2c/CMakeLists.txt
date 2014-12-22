# test: test_2cstokes2p2c_1ff
add_dumux_test(test_2cstokes2p2c_1ff test_2cstokes2p2c test_2cstokes2p2c.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/2cstokes2p2c-ff-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/stokes2c-00007.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_2cstokes2p2c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_2cstokes2p2c_reference.input  
  )

# test: test_2cstokes2p2c_2pm
add_dumux_test(test_2cstokes2p2c_2pm test_2cstokes2p2c test_2cstokes2p2c.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/2cstokes2p2c-pm-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/2p2c-00007.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_2cstokes2p2c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_2cstokes2p2c_reference.input  
  )
