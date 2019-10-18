dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_projection_2d1d
              SOURCES test_projection_2d1d.cc
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-functions_FOUND )"
              LABELS unit discretization)
