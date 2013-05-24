# build target for the CO2 test problem
add_executable("test_boxco2" test_boxco2.cc)
target_link_libraries(test_boxco2 ${DUNE_LIBS})
add_dune_alugrid_flags(test_boxco2)
add_test(test_boxco2
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/co2box-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/heterogeneousbox-00020.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_boxco2
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_boxco2.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/heterogeneousSmall.dgf)

add_executable("test_ccco2" test_ccco2.cc)
add_dune_alugrid_flags(test_ccco2)
add_test(test_ccco2
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/co2cc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/heterogeneouscc-00019.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_ccco2
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_ccco2.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/heterogeneousSmall.dgf)

add_test(test_restartco2
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/co2box-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/heterogeneousbox-00020.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_boxco2
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_boxco2.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/heterogeneousSmall.dgf
  -Restart 42008.6
  -TimeManager.DtInitial 23452.1)

set(TESTS
  test_boxco2
  test_ccco2
  test_restart)

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS})
