add_subdirectory(nonisothermal)

add_input_file_links()

dune_add_test(NAME test_pnm_1p_noncreeping_flow
              SOURCES main.cc
              LABELS porenetworkflow
              COMPILE_DEFINITIONS ISOTHERMAL=1
              COMMAND ${dumux_INCLUDE_DIRS}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_1p_noncreeping_flow-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1p_noncreeping_flow-00000.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1p_noncreeping_flow")
