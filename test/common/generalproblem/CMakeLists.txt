# build target for the simple twophase lens problem in a general
# definition for both box and decoupled model
add_executable("test_generalproblem2p" test_generalproblem2p.cc)

set(NORMALTESTS test_generalproblem2p)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
  add_test(${_exe} ${_exe})
endforeach(_exe ${ALLTESTS})