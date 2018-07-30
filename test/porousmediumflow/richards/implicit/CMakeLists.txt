add_input_file_links()

# isothermal tests
add_executable(test_richards_tpfa EXCLUDE_FROM_ALL test_richardslens_fv.cc)
target_compile_definitions(test_richards_tpfa PUBLIC TYPETAG=RichardsLensCCTypeTag)

add_executable(test_richards_box EXCLUDE_FROM_ALL test_richardslens_fv.cc)
target_compile_definitions(test_richards_box PUBLIC TYPETAG=RichardsLensBoxTypeTag)

dune_add_test(TARGET test_richards_box
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardslensbox-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_richards_box-00007.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_richards_box test_richardslens.input -Problem.Name test_richards_box")

dune_add_test(TARGET test_richards_tpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardslenscc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_richards_tpfa-00007.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_richards_tpfa test_richardslens.input -Problem.Name test_richards_tpfa")

dune_add_test(NAME test_richards_box_parallel_yasp
              TARGET test_richards_box
              CMAKE_GUARD MPI_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"process rank":100}
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardslensbox-reference-parallel.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0000-test_richards_box_parallel_yasp-00007.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/richardslensbox-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_richards_box_parallel_yasp-00007.pvtu
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_richards_box test_richardslens.input -Problem.Name test_richards_box_parallel_yasp -Grid.Overlap 0")

dune_add_test(NAME test_richards_tpfa_parallel_yasp
              TARGET test_richards_tpfa
              CMAKE_GUARD MPI_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"process rank":100}
                      --files ${CMAKE_SOURCE_DIR}/test/references/richardslenscc-reference-parallel.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0000-test_richards_tpfa_parallel_yasp-00007.vtu
                              ${CMAKE_SOURCE_DIR}/test/references/richardslenscc-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_richards_tpfa_parallel_yasp-00007.pvtu
                      --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_richards_tpfa test_richardslens.input -Problem.Name test_richards_tpfa_parallel_yasp -Grid.Overlap 1")

# TODO: there is a small difference because of the precision loss when writing and reading to and from vtk with single precision
# This is why we need a bit higher saturation threshold since pc-sw is very sensitive to pressure for saturations close to 0
dune_add_test(NAME test_richards_tpfa_parallel_yasp_restart
              TARGET test_richards_tpfa
              CMAKE_GUARD MPI_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"process rank":100,"S_w":1e-3,"water content":1e-3}
                      --files ${CMAKE_SOURCE_DIR}/test/references/richardslenscc-reference-parallel.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0000-test_richards_tpfa_parallel_yasp_restart-00004.vtu
                              ${CMAKE_SOURCE_DIR}/test/references/richardslenscc-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_richards_tpfa_parallel_yasp_restart-00004.pvtu
                      --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_richards_tpfa test_richardslens.input -Problem.Name test_richards_tpfa_parallel_yasp_restart -Grid.Overlap 1 -Restart.Time 536.797 -Restart.File s0002-test_richards_tpfa_parallel_yasp-00003.pvtu -TimeLoop.DtInitial 408.68")

# the restart test has to run after the test that produces the corresponding vtu file
set_tests_properties(test_richards_tpfa_parallel_yasp_restart PROPERTIES DEPENDS test_richards_tpfa_parallel_yasp)

dune_add_test(NAME test_richards_box_parallel_ug
              SOURCES test_richardslens_fv.cc
              CMAKE_GUARD "( MPI_FOUND AND HAVE_UG )"
              COMPILE_DEFINITIONS TYPETAG=RichardsLensBoxTypeTag GRIDTYPE=Dune::UGGrid<2>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"process rank":100}
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardslensbox-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_richards_box_parallel_ug-00007.pvtu
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_richards_box_parallel_ug test_richardslens.input -Problem.Name test_richards_box_parallel_ug")

dune_add_test(NAME test_richards_tpfa_parallel_ug
              SOURCES test_richardslens_fv.cc
              CMAKE_GUARD "( MPI_FOUND AND HAVE_UG )"
              COMPILE_DEFINITIONS TYPETAG=RichardsLensCCTypeTag GRIDTYPE=Dune::UGGrid<2>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"process rank":100}
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardslenscc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_richards_tpfa_parallel_ug-00007.pvtu
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_richards_tpfa_parallel_ug test_richardslens.input -Problem.Name test_richards_tpfa_parallel_ug")

dune_add_test(NAME test_richards_box_parallel_alu
              SOURCES test_richardslens_fv.cc
              CMAKE_GUARD "( MPI_FOUND AND dune-alugrid_FOUND )"
              COMPILE_DEFINITIONS TYPETAG=RichardsLensBoxTypeTag GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"process rank":100}
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardslensbox-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_richards_box_parallel_alu-00007.pvtu
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_richards_box_parallel_alu test_richardslens.input -Problem.Name test_richards_box_parallel_alu")

dune_add_test(NAME test_richards_tpfa_parallel_alu
              SOURCES test_richardslens_fv.cc
              CMAKE_GUARD "( MPI_FOUND AND dune-alugrid_FOUND )"
              COMPILE_DEFINITIONS TYPETAG=RichardsLensCCTypeTag GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"process rank":100}
                       --files ${CMAKE_SOURCE_DIR}/test/references/richardslenscc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_richards_tpfa_parallel_alu-00007.pvtu
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_richards_tpfa_parallel_alu test_richardslens.input -Problem.Name test_richards_tpfa_parallel_alu")

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
