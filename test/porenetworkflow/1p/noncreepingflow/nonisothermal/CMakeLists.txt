add_input_file_links()

dune_add_test(NAME test_pnm_1pni_noncreeping_flow
              SOURCES main.cc
              LABELS porenetworkflow
              COMPILE_DEFINITIONS ISOTHERMAL=0
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_1pni_noncreeping_flow-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1pni_noncreeping_flow-00014.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1pni_noncreeping_flow")
