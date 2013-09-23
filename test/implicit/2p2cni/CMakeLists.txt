add_dumux_test(test_box2p2cni test_box2p2cni test_box2p2cni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/waterairbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/waterairbox-00010.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box2p2cni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box2p2cni.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_2p2cni.dgf
  -TimeManager.TEnd 1e4
  -TimeManager.DtInitial 250)

add_dumux_test(test_cc2p2cni test_cc2p2cni test_cc2p2cni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/wateraircc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/wateraircc-00010.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc2p2cni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_cc2p2cni.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_2p2cni.dgf
  -TimeManager.TEnd 1e4
  -TimeManager.DtInitial 250)