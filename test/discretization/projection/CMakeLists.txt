dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_l2_projection_2d1d
               SOURCES test_projection_2d1d.cc
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-functions_FOUND )"
               LABELS unit discretization)

dumux_add_test(NAME test_l2_projection
               SOURCES test_l2_projection.cc
               CMAKE_GUARD dune-functions_FOUND
               LABELS unit discretization)
