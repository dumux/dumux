add_input_file_links(FILES params.input)

dune_add_test(NAME test_3pni_conduction_box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=ThreePNIConductionBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_3pni_conduction_box-reference.vtu
                            ${CMAKE_CURRENT_BINARY_DIR}/test_3pni_conduction_box-00005.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_3pni_conduction_box params.input -Problem.Name test_3pni_conduction_box"
                            --zeroThreshold {"velocity_aq \(m/s\)_1":1e-8})

dune_add_test(NAME test_3pni_conduction_tpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=ThreePNIConductionCCTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_3pni_conduction_tpfa-reference.vtu
                            ${CMAKE_CURRENT_BINARY_DIR}/test_3pni_conduction_tpfa-00005.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_3pni_conduction_tpfa params.input -Problem.Name test_3pni_conduction_tpfa"
                            --zeroThreshold {"velocity_aq \(m/s\)_1":1e-8})
