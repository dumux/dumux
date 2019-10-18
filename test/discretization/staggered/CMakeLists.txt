dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_staggeredfvgeometry
              SOURCES test_staggeredfvgeometry.cc
              LABELS unit discretization)

dumux_add_test(NAME test_staggered_free_flow_geometry
              SOURCES test_staggered_free_flow_geometry.cc
              LABELS unit discretization)
