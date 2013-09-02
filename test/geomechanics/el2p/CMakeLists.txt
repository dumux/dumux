# add build targets
add_executable("test_el2p" EXCLUDE_FROM_ALL test_el2p.cc)
add_dune_superlu_flags(test_el2p)
add_test(test_el2p
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/el2p-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/el2p-00003.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_el2p
  ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_el2p.dgf 100 200 100)

set(NORMALTESTS test_el2p)
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
