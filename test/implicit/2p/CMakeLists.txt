add_dumux_test(test_box2p test_box2p test_box2p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/lensbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/lensbox-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box2p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box2p.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_2p.dgf
  -TimeManager.TEnd 3000
  -TimeManager.DtInitial 250)

add_dumux_test(test_cc2p test_cc2p test_cc2p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/lenscc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/lenscc-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc2p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_cc2p.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_2p.dgf
  -TimeManager.TEnd 3000
  -TimeManager.DtInitial 250)
