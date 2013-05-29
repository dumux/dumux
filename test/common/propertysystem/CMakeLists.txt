# build the test for the property system
add_executable("test_propertysystem" EXCLUDE_FROM_ALL test_propertysystem.cc)
add_test(test_propertysystem
  ${CMAKE_CURRENT_BINARY_DIR}/test_propertysystem)

set(NORMALTESTS test_propertysystem)
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