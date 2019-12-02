add_input_file_links("params.input")
dune_symlink_to_source_files(FILES grids)

dumux_add_test(NAME test_1p_pointsources_timeindependent_tpfa
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              COMPILE_DEFINITIONS TYPETAG=OnePSingularityCCTpfa
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_pointsources_timeindependent_cc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1p_pointsources_timeindependent_tpfa-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_pointsources_timeindependent_tpfa params.input -Problem.Name test_1p_pointsources_timeindependent_tpfa")

dumux_add_test(NAME test_1p_pointsources_timeindependent_box
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              COMPILE_DEFINITIONS TYPETAG=OnePSingularityBox
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_pointsources_timeindependent_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1p_pointsources_timeindependent_box-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_pointsources_timeindependent_box params.input -Problem.Name test_1p_pointsources_timeindependent_box")

dumux_add_test(NAME test_1p_pointsources_timeindependent_tpfa_prism
              LABELS porousmediumflow 1p
              SOURCES main.cc
              CMAKE_GUARD dune-uggrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::UGGrid<3>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              COMPILE_DEFINITIONS TYPETAG=OnePSingularityCCTpfa
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_pointsources_timeindependent_tpfa_prism-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1p_pointsources_timeindependent_tpfa_prism-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_pointsources_timeindependent_tpfa_prism params.input
                                                                            -Problem.Name test_1p_pointsources_timeindependent_tpfa_prism
                                                                            -Grid.File ./grids/well_prism.msh")

dumux_add_test(NAME test_1p_pointsources_timeindependent_box_prism
              LABELS porousmediumflow 1p
              SOURCES main.cc
              CMAKE_GUARD dune-uggrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::UGGrid<3>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              COMPILE_DEFINITIONS TYPETAG=OnePSingularityBox
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_pointsources_timeindependent_box_prism-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1p_pointsources_timeindependent_box_prism-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_pointsources_timeindependent_box_prism params.input
                                                                            -Problem.Name test_1p_pointsources_timeindependent_box_prism
                                                                            -Grid.File ./grids/well_prism.msh")

#install sources
install(FILES
problem.hh
../spatialparams.hh
main.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/1p/implicit/pointsources)
