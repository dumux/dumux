add_subdirectory(reference_data)
dune_symlink_to_source_files(FILES "params_re1.input" "params_re1000.input" "run_and_plot.py")

dumux_add_test(NAME example_ff_liddrivencavity
               SOURCES main.cc
               LABELS freeflow navierstokes example
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_closedsystem_ldc_re1-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/example_ff_liddrivencavity_re1-00002.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/example_ff_liddrivencavity params_re1.input
                             -Grid.Cells \"64 64\"")

dumux_add_test(NAME example_ff_liddrivencavity_re1000
               TARGET example_ff_liddrivencavity
               LABELS freeflow navierstokes example
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_closedsystem_ldc_re1000-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/example_ff_liddrivencavity_re1000-00009.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/example_ff_liddrivencavity params_re1000.input
                             -Grid.Cells \"64 64\" -TimeLoop.TEnd 50")

# test plot script
dune_symlink_to_source_files(run_and_plot.py)
dune_python_add_test(NAME example_ff_liddrivencavity_plot
                     COMMAND ${PYTHON_EXECUTABLE} run_and_plot.py -s -n
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(example_ff_liddrivencavity_plot
                     PROPERTIES LABELS example
                     DEPENDS "example_ff_liddrivencavity;example_ff_liddrivencavity_re1000")

set_tests_properties(example_ff_liddrivencavity_plot PROPERTIES SKIP_RETURN_CODE 127)
