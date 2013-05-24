# add build targets
add_executable("test_impes" test_impes.cc)
add_executable("test_impesadaptive" test_impesadaptive.cc)
add_dune_alugrid_flags(test_impesadaptive)
add_executable("test_impeswithamg" test_impeswithamg.cc)
add_executable("test_transport" test_transport.cc)
add_executable("test_mpfa2p" test_mpfa2p.cc)
add_dune_alugrid_flags(test_mpfa2p)
add_dune_superlu_flags(test_mpfa2p)
>>>>>>> [CMake] All tests are built except the ones in need of PDELab.

set(NORMALTESTS
  test_impes
  test_impesadaptive
  test_transport
  test_mpfa2p)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
  add_test(${_exe} ${_exe})
endforeach(_exe ${ALLTESTS})
