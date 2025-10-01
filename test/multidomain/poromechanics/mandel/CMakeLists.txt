add_subdirectory(analyticalsolution)

dune_symlink_to_source_files(FILES params.input plot_results.py)
dumux_add_test(
    NAME test_md_mandel
    SOURCES test_mandel.cc
    COMMAND ${DUMUX_RUNTEST}
    CMD_ARGS --script fuzzy
             --files ${CMAKE_SOURCE_DIR}/test/references/test_mandel_flow-reference.vtu
                     ${CMAKE_CURRENT_BINARY_DIR}/test_mandel_flow-00005.vtu
                     ${CMAKE_SOURCE_DIR}/test/references/test_mandel_mech-reference.vtu
                     ${CMAKE_CURRENT_BINARY_DIR}/test_mandel_mech-00005.vtu
             --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_mandel params.input -TimeLoop.TEnd 50"
)
