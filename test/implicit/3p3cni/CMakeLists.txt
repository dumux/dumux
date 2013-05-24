add_executable("test_box3p3cni" test_box3p3cni.cc)
add_executable("test_cc3p3cni" test_cc3p3cni.cc)

set(NORMALTESTS
  test_box3p3cni
  test_cc3p3cni)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
  add_test(${_exe} ${_exe})
endforeach(_exe ${ALLTESTS})
