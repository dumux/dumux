add_executable("test_box3p3c" test_box3p3c.cc)
add_test(test_box3p3c
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/infiltrationbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/infiltrationbox-00007.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box3p3c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box3p3c.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_3p3c_coarse.dgf
  -TimeManager.TEnd 8.64e5
  -TimeManager.DtInitial 8.64e4)

add_executable("test_cc3p3c" test_cc3p3c.cc)
add_test(test_cc3p3c
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/infiltrationcc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/infiltrationcc-00006.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc3p3c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_cc3p3c.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_3p3c_coarse.dgf
  -TimeManager.TEnd 8.64e5
  -TimeManager.DtInitial 8.64e4)


set(NORMALTESTS
  test_box3p3c
  test_cc3p3c)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
endforeach(_exe ${ALLTESTS})
