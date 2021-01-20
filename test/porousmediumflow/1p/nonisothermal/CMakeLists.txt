add_input_file_links()

# non-isothermal tests
dumux_add_test(NAME test_1pni_conduction_box
              LABELS porousmediumflow 1p  1pni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNIConductionBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1pni_conduction_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1pni_conduction_box-00005.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pni_conduction_box params_conduction.input -Problem.Name test_1pni_conduction_box"
                       --zeroThreshold {"velocity_liq":1e-8})

add_executable(test_1pni_convection_box EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_1pni_convection_box PUBLIC TYPETAG=OnePNIConvectionBox)

dumux_add_test(NAME test_1pni_convection_box
              TARGET test_1pni_convection_box
              LABELS porousmediumflow 1p  1pni
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1pni_convection_box-reference.vtp
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1pni_convection_box-00009.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pni_convection_box params_convection.input -Problem.Name test_1pni_convection_box"
                       --zeroThreshold {"velocity":1e-15})

dumux_add_test(NAME test_1pni_convection_box_restart
              TARGET test_1pni_convection_box
              LABELS porousmediumflow 1p  1pni
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1pni_convection_box-reference.vtp
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1pni_convection_box_restart-00004.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pni_convection_box params_convection.input -Problem.Name test_1pni_convection_box_restart -TimeLoop.DtInitial 1000 -Restart.Time 16768 -Restart.File test_1pni_convection_box-00006.vtp"
                       --zeroThreshold {"velocity":1e-15})

# the restart test has to run after the test that produces the corresponding vtu file
set_tests_properties(test_1pni_convection_box_restart PROPERTIES DEPENDS test_1pni_convection_box)

dumux_add_test(NAME test_1pni_conduction_tpfa
              LABELS porousmediumflow 1p  1pni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNIConductionCCTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1pni_conduction_cc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1pni_conduction_tpfa-00005.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pni_conduction_tpfa params_conduction.input -Problem.Name test_1pni_conduction_tpfa"
                       --zeroThreshold {"velocity":1e-8})

dumux_add_test(NAME test_1pni_convection_tpfa
              LABELS porousmediumflow 1p  1pni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNIConvectionCCTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1pni_convection_cc-reference.vtp
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1pni_convection_tpfa-00009.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pni_convection_tpfa params_convection.input -Problem.Name test_1pni_convection_tpfa")

dumux_add_test(NAME test_1pni_conduction_mpfa
              LABELS porousmediumflow 1p  1pni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNIConductionCCMpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1pni_conduction_mpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1pni_conduction_mpfa-00005.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pni_conduction_mpfa params_conduction.input -Problem.Name test_1pni_conduction_mpfa")
