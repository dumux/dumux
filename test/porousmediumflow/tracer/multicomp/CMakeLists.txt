add_input_file_links()

# tracer tests
dune_add_test(NAME test_maxwellstefan_tpfa
              SOURCES test_tracer_maxwellstefan.cc
              COMPILE_DEFINITIONS TYPETAG=MaxwellStefanTestCCTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/cctracermaxwellstefan-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_ccmaxwellstefan-00026.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_maxwellstefan_tpfa test_tracer_maxwellstefan.input -Problem.Name test_ccmaxwellstefan")

dune_add_test(NAME test_maxwellstefan_box
              SOURCES test_tracer_maxwellstefan.cc
              COMPILE_DEFINITIONS TYPETAG=MaxwellStefanTestBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/boxtracermaxwellstefan-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_boxmaxwellstefan-00026.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_maxwellstefan_box test_tracer_maxwellstefan.input -Problem.Name test_boxmaxwellstefan")
#install sources
install(FILES
maxwellstefantestproblem.hh
maxwellstefantestspatialparams.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/tracer)
