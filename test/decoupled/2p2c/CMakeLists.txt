add_dumux_test(test_adaptive2p2c test_adaptive2p2c test_adaptive2p2c.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_adaptive2p2c-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_adaptive2p2c-00008.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_adaptive2p2c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_adaptive2p2c.input)
add_dune_alugrid_flags(test_adaptive2p2c)

add_dumux_test(test_dec2p2c test_dec2p2c test_dec2p2c.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_dec2p2c-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_dec2p2c-00021.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_dec2p2c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_dec2p2c.input)

add_dumux_test(test_multiphysics2p2c test_multiphysics2p2c test_multiphysics2p2c.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_multiphysics2p2c-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_multiphysics2p2c-00021.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_multiphysics2p2c)
