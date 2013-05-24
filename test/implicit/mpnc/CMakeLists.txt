# build target for the twophase-twocomponent test problem
add_executable("test_boxmpnc" test_boxmpnc.cc)
add_executable("test_ccmpnc" test_ccmpnc.cc)
target_link_libraries("test_ccmpnc" ${DumuxLinkLibraries})

# build target for the one-phase forchheimer test problem
add_executable("test_forchheimer1p" test_forchheimer1p.cc)

# build target for the two-phase forchheimer test problem
add_executable("test_forchheimer2p" test_forchheimer2p.cc)

set(NORMALTESTS
  test_boxmpnc
  test_ccmpnc
  test_forchheimer1p
  test_forchheimer2p)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
  add_test(${_exe} ${_exe})
endforeach(_exe ${ALLTESTS})
