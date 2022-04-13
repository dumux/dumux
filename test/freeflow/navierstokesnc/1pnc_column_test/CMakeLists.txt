dumux_add_test(NAME test_ff_stokes2c_column_test
              LABELS freeflow navierstokes
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK)

dune_symlink_to_source_files(FILES "params.input")
