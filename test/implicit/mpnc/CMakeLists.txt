add_dumux_test(test_boxmpnc test_boxmpnc test_boxmpnc.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/obstaclebox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/obstaclebox-00010.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_boxmpnc
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_boxmpnc.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/obstacle_24x16.dgf
  -TimeManager.TEnd 1e4
  -TimeManager.DtInitial 250)

add_dumux_test(test_ccmpnc test_ccmpnc test_ccmpnc.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/obstaclecc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/obstaclecc-00010.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_ccmpnc
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_ccmpnc.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/obstacle_24x16.dgf
  -TimeManager.TEnd 1e4
  -TimeManager.DtInitial 250)

# build target for the one-phase Forchheimer test problem
add_dumux_test(test_forchheimer1p test_forchheimer1p test_forchheimer1p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/forchheimer1p-reference.vtp
  ${CMAKE_CURRENT_BINARY_DIR}/test_forchheimer1p-00002.vtp
  ${CMAKE_CURRENT_BINARY_DIR}/test_forchheimer1p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_forchheimer1p.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/forchheimer1d.dgf)

# build target for the two-phase Forchheimer test problem
add_dumux_test(test_forchheimer2p test_forchheimer2p test_forchheimer2p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/forchheimer2p-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_forchheimer2p-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_forchheimer2p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_forchheimer2p.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/obstacle_24x16.dgf)

# build target for the full kinetic test problem
add_dumux_test(test_boxmpnckinetic test_boxmpnckinetic test_boxmpnckinetic.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/evaporationatmosphere-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/evaporationatmosphere-00006.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_boxmpnckinetic
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_boxmpnckinetic.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/interfacedomain.dgf)
  
  # build target for the energy kinetic test problem, two energy balance equations
add_dumux_test(test_boxmpncthermalnonequil test_boxmpncthermalnonequil test_boxmpncthermalnonequil.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/combustion-reference.vtp
  ${CMAKE_CURRENT_BINARY_DIR}/combustion-00078.vtp
  ${CMAKE_CURRENT_BINARY_DIR}/test_boxmpncthermalnonequil
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_boxmpncthermalnonequil.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/combustionOutflowGridLinNX100LogNx100.dgf)
add_dune_superlu_flags(test_boxmpncthermalnonequil)

  
target_link_libraries(test_boxmpnckinetic ${DUNE_LIBS}  ${MPI_LIBRARIES})
add_dune_alugrid_flags(test_boxmpnckinetic)
add_dune_ug_flags(test_boxmpnckinetic)
