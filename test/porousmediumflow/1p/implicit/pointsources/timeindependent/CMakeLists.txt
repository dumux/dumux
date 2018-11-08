add_input_file_links("params.input")

dune_add_test(NAME test_1p_pointsources_timeindependent_tpfa
              SOURCES main.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              COMPILE_DEFINITIONS TYPETAG=OnePSingularityCCTpfa
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_pointsources_timeindependent_cc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1p_pointsources_timeindependent_tpfa-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_pointsources_timeindependent_tpfa params.input -Problem.Name test_1p_pointsources_timeindependent_tpfa")

dune_add_test(NAME test_1p_pointsources_timeindependent_box
              SOURCES main.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              COMPILE_DEFINITIONS TYPETAG=OnePSingularityBox
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_pointsources_timeindependent_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1p_pointsources_timeindependent_box-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_pointsources_timeindependent_box params.input -Problem.Name test_1p_pointsources_timeindependent_box")


#install sources
install(FILES
problem.hh
../spatialparams.hh
main.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/1p/implicit/pointsources)
