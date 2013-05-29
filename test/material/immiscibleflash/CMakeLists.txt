add_executable("test_immiscibleflash" EXCLUDE_FROM_ALL test_immiscibleflash.cc)
add_test(test_immiscibleflash
  ${CMAKE_CURRENT_BINARY_DIR}/test_immiscibleflash)

set(NORMALTESTS test_immiscibleflash)
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
