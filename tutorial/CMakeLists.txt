# add build targets
add_executable("tutorial_decoupled" tutorial_decoupled.cc)
add_executable("tutorial_coupled" tutorial_coupled.cc)
add_dune_alugrid_flags(tutorial_coupled)

set(NORMALTESTS
  tutorial_decoupled
  tutorial_coupled)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
  add_test(${_exe} ${_exe})
endforeach(_exe ${ALLTESTS})
