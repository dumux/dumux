# build target for the 2p2cni test problem
add_executable("test_box2p2cni" test_box2p2cni.cc)
add_executable("test_cc2p2cni" test_cc2p2cni.cc)

set(NORMALTESTS
  test_box2p2cni
  test_cc2p2cni)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
  add_test(${_exe} ${_exe})
endforeach(_exe ${ALLTESTS})
