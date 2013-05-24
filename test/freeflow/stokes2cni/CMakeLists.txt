# add build targets
add_executable("test_stokes2cni" test_stokes2cni.cc)
add_dune_superlu_flags(test_stokes2cni)
add_test(test_stokes
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/stokes2cni-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/stokes2cni-00008.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_stokes2cni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_stokes2cni.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_stokes2cni.dgf)

set(NORMALTESTS test_stokes2cni)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
endforeach(_exe ${ALLTESTS})
