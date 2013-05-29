# build target for the onephase-onecomponent test problem
add_executable("test_box1p" EXCLUDE_FROM_ALL test_box1p.cc)
add_test(test_box1p
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/1ptestbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/1ptestbox-00002.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box1p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box1p.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_1p_2d.dgf               
  -TimeManager.TEnd 1
  -TimeManager.DtInitial 1)

add_executable("test_box1pwithamg" EXCLUDE_FROM_ALL test_box1pwithamg.cc)
add_test(test_box1pwithamg
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/1ptestbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/1ptestboxwithamg-00002.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box1pwithamg
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box1pwithamg.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_1p_2d.dgf
  -TimeManager.TEnd 1
  -TimeManager.DtInitial 1)

add_executable("test_cc1p" EXCLUDE_FROM_ALL test_cc1p.cc)
add_test(test_cc1p
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/1ptestcc-00002.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc1p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_cc1p.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_1p_2d.dgf
  -TimeManager.TEnd 1
  -TimeManager.DtInitial 1)

add_executable("test_cc1pwithamg" EXCLUDE_FROM_ALL test_cc1pwithamg.cc)
add_test(test_cc1pwithamg
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/1ptestccwithamg-00002.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc1pwithamg
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_cc1pwithamg.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_1p_2d.dgf
  -TimeManager.TEnd 1
  -TimeManager.DtInitial 1)

set(NORMALTESTS
  test_box1p
  test_box1pwithamg
  test_cc1p
  test_cc1pwithamg)
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
