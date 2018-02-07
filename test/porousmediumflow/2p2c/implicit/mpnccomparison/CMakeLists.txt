add_input_file_links()

dune_add_test(SOURCES test_2p2c_comparison_fv.cc
              NAME test_2p2c_comparison_box
              COMPILE_DEFINITIONS TYPETAG=TwoPTwoCComparisonBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/mpnc_2p2c_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_comparison_box-00009.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_comparison_box test_2p2c_comparison_fv.input -Problem.Name test_2p2c_comparison_box")

dune_add_test(SOURCES test_2p2c_comparison_fv.cc
              NAME test_2p2c_comparison_tpfa
              COMPILE_DEFINITIONS TYPETAG=TwoPTwoCComparisonCCTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/mpnc_2p2c_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_comparison_tpfa-00009.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_comparison_tpfa test_2p2c_comparison_fv.input -Problem.Name test_2p2c_comparison_tpfa")

#install sources
install(FILES
2p2c_comparison_problem.hh
2p2c_comparison_spatialparams.hh
test_2p2c_comparison_fv.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2p2c/comparison)
