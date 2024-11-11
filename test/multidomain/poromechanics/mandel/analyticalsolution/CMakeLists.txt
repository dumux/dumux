dune_symlink_to_source_files(FILES params.input pos.txt)

dumux_add_test(
    NAME test_md_mandel_analytical_solution
    SOURCES test_analyticalsolution.cc
    COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
    CMD_ARGS --script fuzzyData --delimiter " "
             --files ${CMAKE_SOURCE_DIR}/test/references/test_md_mandel_p-5000_reference.txt
                     ${CMAKE_CURRENT_BINARY_DIR}/p-5000.txt
                     ${CMAKE_SOURCE_DIR}/test/references/test_md_mandel_ux-5000_reference.txt
                     ${CMAKE_CURRENT_BINARY_DIR}/ux-5000.txt
             --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_mandel_analytical_solution params.input"
)
