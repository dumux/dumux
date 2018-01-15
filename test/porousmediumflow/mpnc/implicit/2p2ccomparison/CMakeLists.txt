add_input_file_links()

dune_add_test(NAME test_mpnc_comparison_box
              SOURCES test_mpnc_comparison_fv.cc
              COMPILE_DEFINITIONS TYPETAG=MPNCComparisonBoxProblem
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/mpnc_2p2c_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_comparison_box-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_comparison_box test_mpnc_comparison_fv.input -Problem.Name test_mpnc_comparison_box")

dune_add_test(NAME test_mpnc_comparison_tpfa
              SOURCES test_mpnc_comparison_fv.cc
              COMPILE_DEFINITIONS TYPETAG=MPNCComparisonCCProblem
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/mpnc_2p2c_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_comparison_tpfa-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_comparison_tpfa test_mpnc_comparison_fv.input -Problem.Name test_mpnc_comparison_tpfa")

#install sources
install(FILES
mpnc_comparison_problem.hh
test_mpnc_comparison_fv.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/mpnc/2p2ccomparison)
