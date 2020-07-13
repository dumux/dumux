add_subdirectory(static)

add_input_file_links()

dune_add_test(NAME test_pnm_2p
              SOURCES main.cc
              LABELS porenetworkflow
              COMPILE_DEFINITIONS ISOTHERMAL=1
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p-00171.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p")

dune_add_test(NAME test_pnm_2pni
              SOURCES main.cc
              LABELS porenetworkflow
              COMPILE_DEFINITIONS ISOTHERMAL=0
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2pni-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pni-00190.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pni params_ni.input -Problem.Name test_pnm_2pni")
