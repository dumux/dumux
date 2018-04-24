add_input_file_links()

dune_add_test(NAME test_channel_zeroeq2c
              SOURCES test_channel.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_channel_zeroeq2c.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_channel2c-00031.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_channel_zeroeq2c test_channel2c.input")

dune_add_test(NAME test_channel_zeroeq2cni
              SOURCES test_channel.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_channel_zeroeq2cni.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_channel2cni-00035.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_channel_zeroeq2cni test_channel2cni.input")
target_compile_definitions(test_channel_zeroeq2cni PUBLIC "NONISOTHERMAL=1")

add_executable(test_channel_lowrekepsilon2c EXCLUDE_FROM_ALL test_channel.cc)
target_compile_definitions(test_channel_lowrekepsilon2c PUBLIC "LOWREKEPSILON=1")

dune_add_test(NAME test_channel_lowrekepsilon2cni
              SOURCES test_channel.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_channel_lowrekepsilon2cni.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_channel2cni-00033.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_channel_lowrekepsilon2cni test_channel2cni.input")
target_compile_definitions(test_channel_lowrekepsilon2cni PUBLIC "LOWREKEPSILON=1" "NONISOTHERMAL=1")
