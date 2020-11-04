add_subdirectory(reference_data)

dumux_add_test(NAME example_ff_liddrivencavity
               SOURCES main.cc
               LABELS freeflow navierstokes example
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_closedsystem_ldc_re1-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/example_ff_liddrivencavity_re1-00002.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/example_ff_liddrivencavity params_re1.input
                             -Problem.Name example_ff_liddrivencavity_re1
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
                             -Problem.Name example_ff_liddrivencavity_re1000
                             -Grid.Cells \"64 64\" -TimeLoop.TEnd 50")

# test plot script
set(DUNE_REENABLE_ADD_TEST True)
string(CONCAT LDCLOGFILES
              "example_ff_liddrivencavity_re1_vx.log;example_ff_liddrivencavity_re1_vy.log;"
              "example_ff_liddrivencavity_re1000_vx.log;example_ff_liddrivencavity_re1000_vy.log")
add_test(NAME example_ff_liddrivencavity_plot COMMAND ${CMAKE_CURRENT_BINARY_DIR}/run_and_plot.py -s -n)
set_tests_properties(example_ff_liddrivencavity_plot
                     PROPERTIES LABEL example REQUIRED_FILES "${LDCLOGFILES}")
set(DUNE_REENABLE_ADD_TEST False)

dune_symlink_to_source_files(FILES "params_re1.input" "params_re1000.input" "run_and_plot.py")
