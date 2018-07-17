add_input_file_links()

dune_add_test(NAME test_flatplate_kepsilon2c
              SOURCES test_flatplate.cc
              COMPILE_ONLY CMAKE_GUARD HAVE_UMFPACK)
target_compile_definitions(test_flatplate_kepsilon2c PUBLIC "KEPSILON=1")

add_executable(test_flatplate_kepsilon2cni EXCLUDE_FROM_ALL test_flatplate.cc)
target_compile_definitions(test_flatplate_kepsilon2cni PUBLIC "KEPSILON=1" "NONISOTHERMAL=1")

dune_add_test(NAME test_flatplate_kepsilon2cni_twolayer
              TARGET test_flatplate_kepsilon2cni
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_flatplate_kepsilon2cni.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_flatplate2cni-00033.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_flatplate_kepsilon2cni
                                       test_flatplate2cni.input")

dune_add_test(NAME test_flatplate_kepsilon2cni_wallfunction
              TARGET test_flatplate_kepsilon2cni
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_flatplate_kepsilon2cni_wallfunction.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_flatplate2cni_wallfunction-00024.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_flatplate_kepsilon2cni
                                       test_flatplate2cni_wallfunction.input")

dune_add_test(NAME test_flatplate_komega2c
              SOURCES test_flatplate.cc
              COMPILE_ONLY CMAKE_GUARD HAVE_UMFPACK)
target_compile_definitions(test_flatplate_komega2c PUBLIC "KOMEGA=1")

dune_add_test(NAME test_flatplate_komega2cni
              SOURCES test_flatplate.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_flatplate_komega2cni.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_flatplate2cni-00032.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_flatplate_komega2cni
                                       test_flatplate2cni.input
                                       -Newton.TargetSteps 12")
target_compile_definitions(test_flatplate_komega2cni PUBLIC "KOMEGA=1" "NONISOTHERMAL=1")

dune_add_test(NAME test_flatplate_lowrekepsilon2c
              SOURCES test_flatplate.cc
              COMPILE_ONLY CMAKE_GUARD HAVE_UMFPACK)
target_compile_definitions(test_flatplate_lowrekepsilon2c PUBLIC "LOWREKEPSILON=1")

dune_add_test(NAME test_flatplate_lowrekepsilon2cni
              SOURCES test_flatplate.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_flatplate_lowrekepsilon2cni.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_flatplate2cni-00032.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_flatplate_lowrekepsilon2cni test_flatplate2cni.input")
target_compile_definitions(test_flatplate_lowrekepsilon2cni PUBLIC "LOWREKEPSILON=1" "NONISOTHERMAL=1")

dune_add_test(NAME test_flatplate_oneeq2c
              SOURCES test_flatplate.cc
              COMPILE_ONLY CMAKE_GUARD HAVE_UMFPACK)
target_compile_definitions(test_flatplate_oneeq2c PUBLIC "ONEEQ=1")

dune_add_test(NAME test_flatplate_oneeq2cni
              SOURCES test_flatplate.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_flatplate_oneeq2cni.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_flatplate2cni-00031.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_flatplate_oneeq2cni test_flatplate2cni.input")
target_compile_definitions(test_flatplate_oneeq2cni PUBLIC "ONEEQ=1" "NONISOTHERMAL=1")

dune_add_test(NAME test_flatplate_zeroeq2c
              SOURCES test_flatplate.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_flatplate_zeroeq2c.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_flatplate2c-00031.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_flatplate_zeroeq2c test_flatplate2c.input")

dune_add_test(NAME test_flatplate_zeroeq2cni
              SOURCES test_flatplate.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_flatplate_zeroeq2cni.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_flatplate2cni-00032.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_flatplate_zeroeq2cni test_flatplate2cni.input")
target_compile_definitions(test_flatplate_zeroeq2cni PUBLIC "NONISOTHERMAL=1")
