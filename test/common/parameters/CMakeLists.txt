dumux_add_test(SOURCES test_loggingparametertree.cc
               LABELS unit)
dumux_add_test(SOURCES test_paramtranslate.cc
               LABELS unit)
dune_symlink_to_source_files(FILES "params.input" "params_solver.input")
