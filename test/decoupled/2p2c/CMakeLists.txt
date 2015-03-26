add_dumux_test(test_adaptive2p2c2d test_adaptive2p2c2d test_adaptive2p2c2d.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_adaptive2p2c2d-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_adaptive2p2c2d-00008.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_adaptive2p2c2d)

add_dumux_test(test_adaptive2p2c3d test_adaptive2p2c3d test_adaptive2p2c3d.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_adaptive2p2c3d-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_adaptive2p2c3d-00012.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_adaptive2p2c3d)

add_dumux_test(test_dec2p2c test_dec2p2c test_dec2p2c.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_dec2p2c-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_dec2p2c-00021.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_dec2p2c)

add_dumux_test(test_multiphysics2p2c test_multiphysics2p2c test_multiphysics2p2c.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/test_multiphysics2p2c-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_multiphysics2p2c-00021.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_multiphysics2p2c)
