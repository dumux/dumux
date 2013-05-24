if (MPI_FOUND)
  # tell UG that the model is parallelized
  add_definitions(-DModelP)
endif (MPI_FOUND)

# build target for the Richards test problem
add_executable("test_boxrichards" test_boxrichards.cc)
add_dune_ug_flags(test_boxrichards)

add_executable("test_ccrichards" test_ccrichards.cc)
add_dune_ug_flags(test_ccrichards)

set(NORMALTESTS
  test_boxrichards
  test_ccrichards)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
  add_test(${_exe} ${_exe})
endforeach(_exe ${ALLTESTS})
