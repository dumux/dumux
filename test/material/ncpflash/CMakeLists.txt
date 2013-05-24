add_executable("test_ncpflash" test_ncpflash.cc)
add_test(test_ncpflash
  ${CMAKE_CURRENT_BINARY_DIR}/test_ncpflash)

set(NORMALTESTS test_ncpflash)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
endforeach(_exe ${ALLTESTS})
