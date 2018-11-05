add_input_file_links()

dune_add_test(NAME test_ff_rans2c_kepsilon
              SOURCES main.cc
              COMPILE_ONLY CMAKE_GUARD HAVE_UMFPACK)
target_compile_definitions(test_ff_rans2c_kepsilon PUBLIC "KEPSILON=1")

add_executable(test_ff_rans2cni_kepsilon EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_ff_rans2cni_kepsilon PUBLIC "KEPSILON=1" "NONISOTHERMAL=1")

dune_add_test(NAME test_ff_rans2cni_kepsilon_twolayer
              TARGET test_ff_rans2cni_kepsilon
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans2cni_kepsilon_twolayer-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_kepsilon_twolayer-00033.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_kepsilon params_nonisothermal.input
                             -Problem.Name test_ff_rans2cni_kepsilon_twolayer")

dune_add_test(NAME test_ff_rans2cni_kepsilon_wallfunction
              TARGET test_ff_rans2cni_kepsilon
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans2cni_kepsilon_wallfunction-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_kepsilon_wallfunction-00024.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_kepsilon params_niwallfunction.input
                             -Problem.Name test_ff_rans2cni_kepsilon_wallfunction")

dune_add_test(NAME test_ff_rans2c_komega
              SOURCES main.cc
              COMPILE_ONLY CMAKE_GUARD HAVE_UMFPACK)
target_compile_definitions(test_ff_rans2c_komega PUBLIC "KOMEGA=1")

dune_add_test(NAME test_ff_rans2cni_komega
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans2cni_komega-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_komega-00032.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_komega params_nonisothermal.input
                                       -Problem.Name test_ff_rans2cni_komega
                                       -Newton.TargetSteps 12")
target_compile_definitions(test_ff_rans2cni_komega PUBLIC "KOMEGA=1" "NONISOTHERMAL=1")

dune_add_test(NAME test_ff_rans2c_lowrekepsilon
              SOURCES main.cc
              COMPILE_ONLY CMAKE_GUARD HAVE_UMFPACK)
target_compile_definitions(test_ff_rans2c_lowrekepsilon PUBLIC "LOWREKEPSILON=1")

dune_add_test(NAME test_ff_rans2cni_lowrekepsilon
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans2cni_lowrekepsilon-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_lowrekepsilon-00032.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_lowrekepsilon params_nonisothermal.input
                             -Problem.Name test_ff_rans2cni_lowrekepsilon")
target_compile_definitions(test_ff_rans2cni_lowrekepsilon PUBLIC "LOWREKEPSILON=1" "NONISOTHERMAL=1")

dune_add_test(NAME test_ff_rans2c_oneeq
              SOURCES main.cc
              COMPILE_ONLY CMAKE_GUARD HAVE_UMFPACK)
target_compile_definitions(test_ff_rans2c_oneeq PUBLIC "ONEEQ=1")

dune_add_test(NAME test_ff_rans2cni_oneeq
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans2cni_oneeq-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_oneeq-00031.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_oneeq params_nonisothermal.input
                             -Problem.Name test_ff_rans2cni_oneeq")
target_compile_definitions(test_ff_rans2cni_oneeq PUBLIC "ONEEQ=1" "NONISOTHERMAL=1")

dune_add_test(NAME test_ff_rans2c_zeroeq
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans2c_zeroeq.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2c_zeroeq-00031.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2c_zeroeq params.input
                             -Problem.Name test_ff_rans2c_zeroeq")

dune_add_test(NAME test_ff_rans2cni_zeroeq
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans2cni_zeroeq.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_zeroeq-00032.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans2cni_zeroeq params_nonisothermal.input
                             -Problem.Name test_ff_rans2cni_zeroeq")
target_compile_definitions(test_ff_rans2cni_zeroeq PUBLIC "NONISOTHERMAL=1")
