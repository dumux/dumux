add_input_file_links(FILES params.input)

# non-isothermal tests
dumux_add_test(NAME test_2p2cni_evaporation_box
              LABELS porousmediumflow 2p2c  2p2cni
              SOURCES main.cc
              TIMEOUT 3600
              COMPILE_DEFINITIONS TYPETAG=EvaporationConstantComponentBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p2cni_evaporation_constant_component_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2cni_evaporation_constant_component_box-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2cni_evaporation_box params.input -Problem.Name test_2p2cni_evaporation_constant_component_box")


dumux_add_test(NAME test_2p2cni_evaporation_tpfa
              LABELS porousmediumflow 2p2c  2p2cni
              SOURCES main.cc
              TIMEOUT 3600
              COMPILE_DEFINITIONS TYPETAG=EvaporationConstantComponentCCTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p2cni_evaporation_constant_component_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2cni_evaporation_constant_component_tpfa-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2cni_evaporation_tpfa params.input -Problem.Name test_2p2cni_evaporation_constant_component_tpfa")