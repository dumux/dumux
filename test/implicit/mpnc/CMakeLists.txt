add_dumux_test(test_boxmpnc test_boxmpnc test_boxmpnc.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/obstaclebox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/obstaclebox-00010.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_boxmpnc)

add_dumux_test(test_ccmpnc test_ccmpnc test_ccmpnc.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/obstaclecc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/obstaclecc-00010.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_ccmpnc)

# build target for the one-phase Forchheimer test problem
add_dumux_test(test_forchheimer1p test_forchheimer1p test_forchheimer1p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/forchheimer1p-reference.vtp
  ${CMAKE_CURRENT_BINARY_DIR}/test_forchheimer1p-00002.vtp
  ${CMAKE_CURRENT_BINARY_DIR}/test_forchheimer1p)

# build target for the two-phase Forchheimer test problem
add_dumux_test(test_forchheimer2p test_forchheimer2p test_forchheimer2p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/forchheimer2p-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_forchheimer2p-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_forchheimer2p)

# build target for the full kinetic test problem
add_dumux_test(test_boxmpnckinetic test_boxmpnckinetic test_boxmpnckinetic.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/evaporationatmosphere-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/evaporationatmosphere-00006.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_boxmpnckinetic)

  # build target for the energy kinetic test problem, two energy balance equations
add_dumux_test(test_boxmpncthermalnonequil test_boxmpncthermalnonequil test_boxmpncthermalnonequil.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/combustion-reference.vtp
  ${CMAKE_CURRENT_BINARY_DIR}/combustion-00078.vtp
  ${CMAKE_CURRENT_BINARY_DIR}/test_boxmpncthermalnonequil)
