add_input_file_links()

add_executable(test_richards_tpfa EXCLUDE_FROM_ALL test_richardslens_fv.cc)
target_compile_definitions(test_richards_tpfa PUBLIC TYPETAG=RichardsLensCCTypeTag)

add_executable(test_richards_box EXCLUDE_FROM_ALL test_richardslens_fv.cc)
target_compile_definitions(test_richards_box PUBLIC TYPETAG=RichardsLensBoxTypeTag)

dune_add_test(NAME test_richards_box_parallel
              TARGET test_richards_box
              CMAKE_GUARD MPI_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"process rank":100}
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardslensbox-reference-parallel.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0000-test_richards_box_parallel-00008.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/richardslensbox-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_richards_box_parallel-00008.pvtu
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_richards_box test_richardslens.input -Problem.Name test_richards_box_parallel -Grid.Overlap 0")


dune_add_test(NAME test_richards_tpfa_parallel
              TARGET test_richards_tpfa
              CMAKE_GUARD MPI_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"process rank":100}
                      --files ${CMAKE_SOURCE_DIR}/test/references/richardslenscc-reference-parallel.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0000-test_richards_tpfa_parallel-00007.vtu
                              ${CMAKE_SOURCE_DIR}/test/references/richardslenscc-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_richards_tpfa_parallel-00007.pvtu
                      --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_richards_tpfa test_richardslens.input -Problem.Name test_richards_tpfa_parallel -Grid.Overlap 1")

# isothermal tests
dune_add_test(TARGET test_richards_box
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardslensbox-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_richards_box-00008.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_richards_box test_richardslens.input -Problem.Name test_richards_box")

dune_add_test(TARGET test_richards_tpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardslenscc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_richards_tpfa-00007.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_richards_tpfa test_richardslens.input -Problem.Name test_richards_tpfa")

# comparison to analytical solution - only with cc
dune_add_test(SOURCES test_ccrichardsanalytical.cc
              NAME test_ccrichardsanalytical
              COMPILE_DEFINITIONS TYPETAG=RichardsAnalyticalCCTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardsanalyticalcc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsanalytical-00002.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsanalytical test_ccrichardsanalytical.input -Problem.Name test_ccrichardsanalytical")
# non-isothermal tests
dune_add_test(SOURCES test_richardsniconvection_fv.cc
              NAME test_boxrichardsniconvection
              COMPILE_DEFINITIONS TYPETAG=RichardsNIConvectionBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardsniconvectionbox-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsniconvection-00044.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsniconvection test_richardsniconvection.input -Problem.Name test_boxrichardsniconvection"
                      )

dune_add_test(SOURCES test_richardsniconvection_fv.cc
              NAME test_ccrichardsniconvection
              COMPILE_DEFINITIONS TYPETAG=RichardsNIConvectionCCTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardsniconvectioncc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsniconvection-00043.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsniconvection test_richardsniconvection.input -Problem.Name test_ccrichardsniconvection")

dune_add_test(SOURCES test_richardsniconduction_fv.cc
              NAME test_boxrichardsniconduction
              COMPILE_DEFINITIONS TYPETAG=RichardsNIConductionBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardsniconductionbox-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsniconduction-00024.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsniconduction test_richardsniconduction.input -Problem.Name test_boxrichardsniconduction"
                        )

dune_add_test(SOURCES test_richardsniconduction_fv.cc
              NAME test_ccrichardsniconduction
              COMPILE_DEFINITIONS TYPETAG=RichardsNIConductionCCTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardsniconductioncc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsniconduction-00024.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsniconduction test_richardsniconduction.input -Problem.Name test_ccrichardsniconduction"
                      )

dune_add_test(SOURCES test_richardsnievaporation_fv.cc
              NAME test_richardsnievaporation_tpfa
              COMPILE_DEFINITIONS TYPETAG=RichardsNIEvaporationCCTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_ccrichardsevaporation-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsnievaporation-00043.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_richardsnievaporation_tpfa test_richardsnievaporation.input -Problem.Name test_ccrichardsnievaporation"
                      )

dune_add_test(SOURCES test_richardsnievaporation_fv.cc
              NAME test_richardsnievaporation_box
              COMPILE_DEFINITIONS TYPETAG=RichardsNIEvaporationBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_boxrichardsevaporation-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsnievaporation-00043.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_richardsnievaporation_box test_richardsnievaporation.input -Problem.Name test_boxrichardsnievaporation"
                      )


#install sources
install(FILES
richardsanalyticalproblem.hh
richardsanalyticalspatialparams.hh
richardslensproblem.hh
richardslensspatialparams.hh
richardsniconductionproblem.hh
richardsniconvectionproblem.hh
richardsnievaporationproblem.hh
richardsnispatialparams.hh
test_richardslens_fv.cc
test_richardsniconduction_fv.cc
test_richardsniconvection_fv.cc
test_richardsnievaporation_fv.cc
test_ccrichardsanalytical.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/richards)
