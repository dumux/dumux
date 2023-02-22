dune_symlink_to_source_files(FILES "params.input")
dumux_add_test(NAME test_linearsolver
               SOURCES test_linearsolver.cc
               LABELS linear unit)
dumux_add_test(NAME test_scalarproducts
               SOURCES test_scalarproducts.cc
               LABELS linear unit)
dumux_add_test(NAME test_dunevectors
               SOURCES test_dunevectors.cc
               LABELS linear unit)
