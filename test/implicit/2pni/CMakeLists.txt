add_executable("test_box2pni" EXCLUDE_FROM_ALL test_box2pni.cc)
add_test(test_box2pni
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/injection2pnibox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/injection2pnibox-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box2pni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box2pni.input
  -TimeManager.TEnd 1e4
  -TimeManager.DtInitial 250)

add_executable("test_cc2pni" EXCLUDE_FROM_ALL test_cc2pni.cc)
add_test(test_cc2pni
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/injection2pnicc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/injection2pnicc-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc2pni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_cc2pni.input
  -TimeManager.TEnd 1e4
  -TimeManager.DtInitial 250)

set(NORMALTESTS
  test_box2pni
  test_cc2pni)
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
