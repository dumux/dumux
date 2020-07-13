add_input_file_links()

add_executable(test_pnm_1p2c EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_pnm_1p2c PRIVATE ISOTHERMAL=1)

dune_add_test(NAME test_pnm_1p2c_advection
              TARGET test_pnm_1p2c
              LABELS porenetworkflow
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_1p2c-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1p2c-00055.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1p2c")

dune_add_test(NAME test_pnm_1p2c_diffusion
              TARGET test_pnm_1p2c
              LABELS porenetworkflow
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_1p2c_diffusion-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1p2c_diffusion-00054.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1p2c params_diffusion.input")

dune_add_test(NAME test_pnm_1p2cni
              SOURCES main.cc
              LABELS porenetworkflow
              COMPILE_DEFINITIONS ISOTHERMAL=0
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_1p2cni-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1p2cni-00056.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1p2cni -Problem.Name test_pnm_1p2cni")
