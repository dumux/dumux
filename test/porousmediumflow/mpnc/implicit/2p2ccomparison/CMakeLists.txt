add_input_file_links()

dune_add_test(NAME test_mpnc_2p2ccomparison_box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=MPNCComparisonBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_mpnc_2p2ccomparison_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_2p2ccomparison_box-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_2p2ccomparison_box params.input -Problem.Name test_mpnc_2p2ccomparison_box")

dune_add_test(NAME test_mpnc_2p2ccomparison_tpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=MPNCComparisonCC
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_mpnc_2p2ccomparison_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_2p2ccomparison_tpfa-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_2p2ccomparison_tpfa params.input -Problem.Name test_mpnc_2p2ccomparison_tpfa")

#install sources
install(FILES
problem.hh
spatialparams.hh
main.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/mpnc/2p2ccomparison)
