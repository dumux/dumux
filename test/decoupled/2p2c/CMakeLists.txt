# add build targets
add_executable("test_adaptive2p2c" test_adaptive2p2c.cc)
add_dune_alugrid_flags(test_adaptive2p2c)
add_executable("test_dec2p2c" test_dec2p2c.cc)
add_executable("test_multiphysics2p2c" test_multiphysics2p2c.cc)


set(NORMALTESTS
  test_adaptive2p2c
  test_dec2p2c
  test_multiphysics2p2c)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
  add_test(${_exe} ${_exe})
endforeach(_exe ${ALLTESTS})
