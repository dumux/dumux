# build target for the CO2 non-isothermal test problem
add_executable("test_boxco2ni" test_boxco2ni.cc)
add_dune_alugrid_flags(test_boxco2ni)
add_executable("test_ccco2ni" test_ccco2ni.cc)
add_dune_alugrid_flags(test_ccco2ni)

set(NORMALTESTS
  test_boxco2ni
  test_ccco2ni)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
  add_test(${_exe} ${_exe})
endforeach(_exe ${ALLTESTS})
