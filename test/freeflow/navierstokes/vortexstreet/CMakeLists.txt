dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_ff_karman_vortex_street
               SOURCES main.cc
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-subgrid_FOUND )")
