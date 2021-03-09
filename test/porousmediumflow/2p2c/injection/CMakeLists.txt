add_input_file_links(FILES params.input)

dumux_add_test(NAME test_2p2c_injection_box
              LABELS porousmediumflow 2p2c
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=InjectionBox ENABLECACHING=0
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p2c_injection_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_injection_box-00008.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_injection_box params.input -Problem.Name test_2p2c_injection_box")

dumux_add_test(NAME test_2p2c_injection_box_restart
              TARGET test_2p2c_injection_box
              LABELS porousmediumflow 2p2c
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p2c_injection_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_injection_box_restart-00004.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_injection_box params.input -Problem.Name test_2p2c_injection_box_restart -TimeLoop.DtInitial 1488.5 -Restart.Time 2158.85 -Restart.File test_2p2c_injection_box-00004.vtu")

# the restart test has to run after the test that produces the corresponding vtu file
set_tests_properties(test_2p2c_injection_box_restart PROPERTIES DEPENDS test_2p2c_injection_box)

dumux_add_test(NAME test_2p2c_injection_tpfa
              LABELS porousmediumflow 2p2c
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=InjectionCCTpfa ENABLECACHING=0
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p2c_injection_cc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_injection_tpfa-00008.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_injection_tpfa params.input -Problem.Name test_2p2c_injection_tpfa")

dumux_add_test(NAME test_2p2c_injection_mpfa
              LABELS porousmediumflow 2p2c
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=InjectionCCMpfa ENABLECACHING=0
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p2c_injection_cc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_injection_mpfa-00008.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_injection_mpfa params.input -Problem.Name test_2p2c_injection_mpfa")

# isothermal tests with caching
dumux_add_test(NAME test_2p2c_injection_box_caching
              LABELS porousmediumflow 2p2c
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=InjectionBox ENABLECACHING=1
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p2c_injection_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_injection_box_caching-00008.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_injection_box_caching params.input -Problem.Name test_2p2c_injection_box_caching")

dumux_add_test(NAME test_2p2c_injection_tpfa_caching
              LABELS porousmediumflow 2p2c
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=InjectionCCTpfa ENABLECACHING=1
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p2c_injection_cc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_injection_tpfa_caching-00008.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_injection_tpfa_caching params.input -Problem.Name test_2p2c_injection_tpfa_caching")

dumux_add_test(NAME test_2p2c_injection_mpfa_caching
              LABELS porousmediumflow 2p2c
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=InjectionCCMpfa ENABLECACHING=1
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p2c_injection_cc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_injection_mpfa_caching-00008.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_injection_mpfa_caching params.input -Problem.Name test_2p2c_injection_mpfa_caching")
