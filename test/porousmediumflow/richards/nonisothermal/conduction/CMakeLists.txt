add_input_file_links(FILES params.input)

# non-isothermal tests
dumux_add_test(NAME test_richardsni_conduction_box
              LABELS porousmediumflow richards richardsni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=RichardsNIConductionBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_richardsni_conduction_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_richardsni_conduction_box-00024.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_richardsni_conduction_box params.input -Problem.Name test_richardsni_conduction_box")

dumux_add_test(NAME test_richardsni_conduction_tpfa
              LABELS porousmediumflow richards richardsni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=RichardsNIConductionCC
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_richardsni_conduction_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_richardsni_conduction_tpfa-00024.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_richardsni_conduction_tpfa params.input -Problem.Name test_richardsni_conduction_tpfa")
