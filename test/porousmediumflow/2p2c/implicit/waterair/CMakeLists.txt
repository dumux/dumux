add_input_file_links(FILES params.input)

# non-isothermal tests
dune_add_test(NAME test_2p2cni_waterair_box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=WaterAirBox ENABLECACHING=0
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p2cni_waterair_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2cni_waterair_box-00013.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2cni_waterair_box params.input -Problem.Name test_2p2cni_waterair_box")

dune_add_test(NAME test_2p2cni_waterair_buoyancy_box
              TARGET test_2p2cni_waterair_box
              COMPILE_DEFINITIONS TYPETAG=WaterAirBox ENABLECACHING=0
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p2cni_waterair_buoyancy_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2cni_waterair_buoyancy_box-00015.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2cni_waterair_box params.input -Problem.Name test_2p2cni_waterair_buoyancy_box -TimeLoop.DtInitial 150")

dune_add_test(NAME test_2p2cni_waterair_tpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=WaterAirCCTpfa ENABLECACHING=0
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p2cni_waterair_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2cni_waterair_tpfa-00013.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2cni_waterair_tpfa params.input -Problem.Name test_2p2cni_waterair_tpfa")

#install sources
install(FILES
problem.hh
spatialparams.hh
main.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2p2c/waterair)
