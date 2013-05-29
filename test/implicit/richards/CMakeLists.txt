# build target for the Richards test problem
add_executable("test_boxrichards" EXCLUDE_FROM_ALL test_boxrichards.cc)
add_dune_ug_flags(test_boxrichards)
if(MPI_CXX_FOUND)
    # tell UG that the model is parallelized
  add_definitions(-DModelP)
  add_test(test_boxrichards
    ${CMAKE_SOURCE_DIR}/bin/runTest.sh
    ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
    ${CMAKE_SOURCE_DIR}/test/references/richardslensbox-reference-parallel.vtu
    ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0000-richardslensbox-00008.vtu
    "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_boxrichards"
    -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_boxrichards.input
    -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/richardslens-24x16.dgf
    -TimeManager.TEnd 3000
    -TimeManager.DtInitial 100)
else()
  add_test(test_boxrichards
    ${CMAKE_SOURCE_DIR}/bin/runTest.sh
    ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
    ${CMAKE_SOURCE_DIR}/test/references/richardslensbox-reference.vtu
    ${CMAKE_CURRENT_BINARY_DIR}/richardslensbox-00008.vtu
    ${CMAKE_CURRENT_BINARY_DIR}/test_boxrichards
    -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_boxrichards.input
    -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/richardslens-24x16.dgf
    -TimeManager.TEnd 3000
    -TimeManager.DtInitial 100)
endif(MPI_CXX_FOUND)

add_executable("test_ccrichards" EXCLUDE_FROM_ALL test_ccrichards.cc)
add_dune_ug_flags(test_ccrichards)
add_test(test_ccrichards
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/richardslenscc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/richardslenscc-00008.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_ccrichards
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_ccrichards.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/richardslens-24x16.dgf
  -TimeManager.TEnd 3000
  -TimeManager.DtInitial 100)

set(NORMALTESTS
  test_boxrichards
  test_ccrichards)
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
