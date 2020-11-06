add_subdirectory(reference_data)
dune_symlink_to_source_files(FILES "params_re1.input" "params_re1000.input" "run_and_plot.py")
set(DUMUX_EXAMPLE_LDC_CMAKE_GUARDS "( HAVE_UMFPACK )")

dumux_add_test(NAME example_ff_liddrivencavity
               SOURCES main.cc
               LABELS freeflow navierstokes example
               CMAKE_GUARD DUMUX_EXAMPLE_LDC_CMAKE_GUARDS
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_closedsystem_ldc_re1-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/example_ff_liddrivencavity_re1-00002.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/example_ff_liddrivencavity params_re1.input
                             -Grid.Cells \"64 64\"")

dumux_add_test(NAME example_ff_liddrivencavity_re1000
               TARGET example_ff_liddrivencavity
               LABELS freeflow navierstokes example
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_closedsystem_ldc_re1000-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/example_ff_liddrivencavity_re1000-00009.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/example_ff_liddrivencavity params_re1000.input
                             -Grid.Cells \"64 64\" -TimeLoop.TEnd 50")

# test plot script (only if UMFPack is available and the tests above were actually run)
if(DUMUX_EXAMPLE_LDC_CMAKE_GUARDS)
    dumux_add_test(NAME example_ff_liddrivencavity_plot
                   TARGET example_ff_liddrivencavity
                   COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_and_plot.py
                   CMD_ARGS -s -n)

    set_tests_properties(example_ff_liddrivencavity_plot
                         PROPERTIES LABELS example
                         DEPENDS "example_ff_liddrivencavity;example_ff_liddrivencavity_re1000")
endif()
