# build target for the onephase-onecomponent test problem
add_executable("test_box1p" test_box1p.cc)
add_executable("test_box1pwithamg" test_box1pwithamg.cc)
add_executable("test_cc1p" test_cc1p.cc)
add_executable("test_cc1pwithamg" test_cc1pwithamg.cc)

set(NORMALTESTS
  test_box1p
  test_box1pwithamg
  test_cc1p
  test_cc1pwithamg)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
  add_test(${_exe} ${_exe})
endforeach(_exe ${ALLTESTS})
