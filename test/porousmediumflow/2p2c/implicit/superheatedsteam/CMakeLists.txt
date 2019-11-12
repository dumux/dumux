add_input_file_links(FILES params.input)

# non-isothermal tests
dune_add_test(NAME test_2p2cni_superheatedsteam_box
              LABELS porousmediumflow 2p2c  2p2cni
              SOURCES main.cc
              TIMEOUT 1500
              COMPILE_DEFINITIONS TYPETAG=WaterAirBox ENABLECACHING=0
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p2cni_superheatedsteam_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2cni_superheatedsteam_box-00022.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2cni_superheatedsteam_box params.input -Problem.Name test_2p2cni_superheatedsteam_box")

dune_add_test(NAME test_2p2cni_superheatedsteam_tpfa
              LABELS porousmediumflow 2p2c  2p2cni
              SOURCES main.cc
              TIMEOUT 1500
              COMPILE_DEFINITIONS TYPETAG=WaterAirCCTpfa ENABLECACHING=0
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p2cni_superheatedsteam_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2cni_superheatedsteam_tpfa-00015.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2cni_superheatedsteam_tpfa params.input -Problem.Name test_2p2cni_superheatedsteam_tpfa")
