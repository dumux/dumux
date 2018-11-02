add_input_file_links()

dune_add_test(NAME test_1pcctpfa_pointsources
              SOURCES test_1pfv_pointsources.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              COMPILE_DEFINITIONS TYPETAG=OnePSingularityCCTpfa
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/1psingularitycc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/1pcctpfa_singularity-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pcctpfa_pointsources test_1pfv_pointsources.input -Problem.Name 1pcctpfa_singularity")

dune_add_test(NAME test_1pcctpfa_pointsources_timedependent
              SOURCES test_1pfv_pointsources_timedependent.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              COMPILE_DEFINITIONS TYPETAG=OnePSingularityTimeDependentCCTpfa
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/1psingularitycctimedependent-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/1pcctpfa_singularity_timedependent-00003.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pcctpfa_pointsources_timedependent test_1pfv_pointsources.input -TimeLoop.TEnd 4 -Problem.Name 1pcctpfa_singularity_timedependent")

dune_add_test(NAME test_1pbox_pointsources
              SOURCES test_1pfv_pointsources.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              COMPILE_DEFINITIONS TYPETAG=OnePSingularityBox
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/1psingularitybox-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/1pbox_singularity-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pbox_pointsources test_1pfv_pointsources.input -Problem.Name 1pbox_singularity")


#install sources
install(FILES
1psingularityproblem.hh
1psingularityproblemtimedependent.hh
1psingularityspatialparams.hh
test_1p_pointsources.cc
test_1p_pointsources_timedependent.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/1p/implicit/pointsources)
