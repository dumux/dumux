# build target for the CO2 non-isothermal test problem
add_executable("test_boxco2ni" EXCLUDE_FROM_ALL test_boxco2ni.cc)
add_dune_alugrid_flags(test_boxco2ni)
add_test(test_boxco2ni
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/co2nibox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/heterogeneousboxni-00020.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_boxco2ni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_boxco2ni.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/heterogeneousSmall.dgf)

add_executable("test_ccco2ni" EXCLUDE_FROM_ALL test_ccco2ni.cc)
add_dune_alugrid_flags(test_ccco2ni)
add_test(test_ccco2ni
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/co2nicc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/heterogeneousccni-00019.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_ccco2ni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_ccco2ni.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/heterogeneousSmall.dgf)

set(NORMALTESTS
  test_boxco2ni
  test_ccco2ni)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  PROPERTY TEST_INCLUDE_FILE ${CMAKE_CURRENT_BINARY_DIR}/BuildTests.cmake)

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
endforeach(_exe ${ALLTESTS})
