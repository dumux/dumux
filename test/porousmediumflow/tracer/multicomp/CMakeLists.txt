add_input_file_links()

# tracer tests
dune_add_test(NAME test_tracer_maxwellstefan_tpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=MaxwellStefanTestCC
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_tracer_maxwellstefan_tpfa-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_maxwellstefan_tpfa-00026.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_maxwellstefan_tpfa params.input -Problem.Name test_tracer_maxwellstefan_tpfa")

dune_add_test(NAME test_tracer_maxwellstefan_box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=MaxwellStefanTestBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_tracer_maxwellstefan_box-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_maxwellstefan_box-00026.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_maxwellstefan_box params.input -Problem.Name test_tracer_maxwellstefan_box")
#install sources
install(FILES
problem.hh
spatialparams.hh
main.cc
params.input
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/tracer/multicomp)
