dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_facecentered_staggered
              SOURCES test_staggeredfvgeometry.cc
              LABELS unit discretization)