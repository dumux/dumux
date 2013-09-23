# build target for the CO2 non-isothermal test problem
add_dumux_test(test_boxco2ni test_boxco2ni test_boxco2ni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/co2nibox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/heterogeneousboxni-00020.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_boxco2ni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_boxco2ni.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/heterogeneousSmall.dgf)
add_dune_alugrid_flags(test_boxco2ni)

add_dumux_test(test_ccco2ni test_ccco2ni test_ccco2ni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/co2nicc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/heterogeneousccni-00019.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_ccco2ni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_ccco2ni.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/heterogeneousSmall.dgf)
add_dune_alugrid_flags(test_ccco2ni)
