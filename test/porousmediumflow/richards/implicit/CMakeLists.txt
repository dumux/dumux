add_input_file_links()

if(MPI_FOUND)
dune_add_test(SOURCES test_richardslens_fv.cc
              NAME test_boxrichards_parallel
              COMPILE_DEFINITIONS TYPETAG=RichardsLensBoxProblem
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardslensbox-reference-parallel.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0000-richardslensbox-00008.vtu
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_boxrichards test_richardslens.input -Problem.Name test_boxrichards")
endif()
# isothermal tests
dune_add_test(SOURCES test_richardslens_fv.cc
              NAME test_boxrichards
              COMPILE_DEFINITIONS TYPETAG=RichardsLensBoxProblem
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardslensbox-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_boxrichards-00008.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxrichards test_richardslens.input -Problem.Name test_boxrichards")

dune_add_test(SOURCES test_richardslens_fv.cc
              NAME test_ccrichards
              COMPILE_DEFINITIONS TYPETAG=RichardsLensCCProblem
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardslenscc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_ccrichards-00007.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccrichards test_richardslens.input -Problem.Name test_ccrichards")

# comparison to analytical solution - only with cc
dune_add_test(SOURCES test_ccrichardsanalytical.cc
              NAME test_ccrichardsanalytical
              COMPILE_DEFINITIONS TYPETAG=RichardsAnalyticalCCProblem
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardsanalyticalcc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/richardsanalyticalcc-00000.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsanalytical test_ccrichardsanalytical.input -Problem.Name test_ccrichardsanalytical")
# non-isothermal tests
dune_add_test(SOURCES test_richardsniconvection_fv.cc
              NAME test_boxrichardsniconvection
              COMPILE_DEFINITIONS TYPETAG=RichardsNIConvectionBoxProblem
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardsniconvectionbox-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsniconvection-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsniconvection test_richardsniconvection.input -Problem.Name test_boxrichardsniconvection"
                       --zeroThreshold {"velocity":1e-16})

dune_add_test(SOURCES test_richardsniconvection_fv.cc
              NAME test_ccrichardsniconvection
              COMPILE_DEFINITIONS TYPETAG=RichardsNIConvectionCCProblem
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardsniconvectioncc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsniconvection-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsniconvection test_richardsniconvection.input -Problem.Name test_ccrichardsniconvection")

dune_add_test(SOURCES test_richardsniconduction_fv.cc
              NAME test_boxrichardsniconduction
              COMPILE_DEFINITIONS TYPETAG=RichardsNIConductionBoxProblem
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardsniconductionbox-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsniconduction-00006.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxrichardsniconduction test_richardsniconduction.input -Problem.Name test_boxrichardsniconduction"
                       --zeroThreshold {"velocity":1e-8})

dune_add_test(SOURCES test_richardsniconduction_fv.cc
              NAME test_ccrichardsniconduction
              COMPILE_DEFINITIONS TYPETAG=RichardsNIConductionCCProblem
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardsniconductioncc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsniconduction-00006.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccrichardsniconduction test_richardsniconduction.input -Problem.Name test_ccrichardsniconduction"
                       --zeroThreshold {"velocity":1e-8})

#install sources
install(FILES
richardsanalyticalproblem.hh
richardsanalyticalspatialparams.hh
richardslensproblem.hh
richardslensspatialparams.hh
richardsniconductionproblem.hh
richardsniconvectionproblem.hh
richardsnispatialparams.hh
test_richardslens_fv.cc
test_richardsniconduction_fv.cc
test_richardsniconvection_fv.cc
test_ccrichardsanalytical.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/richards)
