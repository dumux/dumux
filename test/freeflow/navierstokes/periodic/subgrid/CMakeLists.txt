add_input_file_links()
dune_symlink_to_source_files(FILES voronoi.pbm)

dumux_add_test(NAME test_periodic_porescale
               SOURCES main.cc
               LABELS porescale periodic freeflow
               CMAKE_GUARD HAVE_UMFPACK)
