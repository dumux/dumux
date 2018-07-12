add_subdirectory(pointsources)
add_subdirectory(incompressible)
add_subdirectory(compressible)
add_subdirectory(periodicbc)

add_input_file_links()
add_gstat_file_links()

dune_symlink_to_source_files(FILES grids)
dune_symlink_to_source_files(FILES "tubesconvergencetest.py")

# isothermal tests
dune_add_test(NAME test_1pcctpfa
              SOURCES test_1pfv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTestCCTpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1ptestcctpfa-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pcctpfa test_1pfv.input -Problem.Name 1ptestcctpfa")

dune_add_test(NAME test_1pccmpfa
              SOURCES test_1pfv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTestCCMpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1ptestccmpfa-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pccmpfa test_1pfv.input -Problem.Name 1ptestccmpfa")

dune_add_test(NAME test_1pbox
              SOURCES test_1pfv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTestBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1ptestbox-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1ptestbox-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pbox test_1pfv.input -Problem.Name 1ptestbox")

# a gstat test (becaue it's a random permeability field we can't test against a reference solution)
dune_add_test(NAME test_1pwithgstat
              SOURCES test_1pfv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTestCCTpfaTypeTag
              CMAKE_GUARD HAVE_GSTAT)

# non-isothermal tests
dune_add_test(NAME test_1pnibox_conduction
              SOURCES test_1pnifv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNIConductionBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/1pniboxconduction-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/1pnibox_conduction-00005.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pnibox_conduction test_1pnifv_conduction.input -Problem.Name 1pnibox_conduction"
                       --zeroThreshold {"velocity_liq":1e-8})

add_executable(test_1pnibox_convection EXCLUDE_FROM_ALL test_1pnifv.cc)
target_compile_definitions(test_1pnibox_convection PUBLIC TYPETAG=OnePNIConvectionBoxTypeTag)

dune_add_test(NAME test_1pnibox_convection
              TARGET test_1pnibox_convection
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/1pniboxconvection-reference.vtp
                               ${CMAKE_CURRENT_BINARY_DIR}/1pnibox_convection-00009.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pnibox_convection test_1pnifv_convection.input -Problem.Name 1pnibox_convection"
                       --zeroThreshold {"velocity":1e-15})

dune_add_test(NAME test_1pnibox_convection_restart
              TARGET test_1pnibox_convection
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/1pniboxconvection-reference.vtp
                               ${CMAKE_CURRENT_BINARY_DIR}/1pnibox_convection_restart-00004.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pnibox_convection test_1pnifv_convection.input -Problem.Name 1pnibox_convection_restart -TimeLoop.DtInitial 1000 -Restart.Time 16768 -Restart.File 1pnibox_convection-00006.vtp"
                       --zeroThreshold {"velocity":1e-15})

# the restart test has to run after the test that produces the corresponding vtu file
set_tests_properties(test_1pnibox_convection_restart PROPERTIES DEPENDS test_1pnibox_convection)

dune_add_test(NAME test_1pnicctpfa_conduction
              SOURCES test_1pnifv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNIConductionCCTpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/1pniccconduction-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/1pnicctpfa_conduction-00005.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pnicctpfa_conduction test_1pnifv_conduction.input -Problem.Name 1pnicctpfa_conduction"
                       --zeroThreshold {"velocity":1e-8})

dune_add_test(NAME test_1pnicctpfa_convection
              SOURCES test_1pnifv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNIConvectionCCTpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/1pniccconvection-reference.vtp
                               ${CMAKE_CURRENT_BINARY_DIR}/1pnicctpfa_convection-00009.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pnicctpfa_convection test_1pnifv_convection.input -Problem.Name 1pnicctpfa_convection")

dune_add_test(NAME test_1pniccmpfa_conduction
              SOURCES test_1pnifv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNIConductionCCMpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/ccmpfa1pniconduction-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/1pniccmpfa_conduction-00005.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pniccmpfa_conduction test_1pnifv_conduction.input -Problem.Name 1pniccmpfa_conduction")

# dim < dimWorld tests with Dune::Foamgrid<1,3>
dune_add_test(NAME test_1pcctpfa_network1d3d
              SOURCES test_1pfv_network1d3d.cc
              COMPILE_DEFINITIONS TYPETAG=TubesTestCCTpfaTypeTag
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ./tubesconvergencetest.py
              CMD_ARGS test_1pcctpfa_network1d3d test_1pfv_network1d3d.input -Problem.Name test_1pcctpfa_network1d3d)

dune_add_test(NAME test_1pbox_network1d3d
              SOURCES test_1pfv_network1d3d.cc
              COMPILE_DEFINITIONS TYPETAG=TubesTestBoxTypeTag
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ./tubesconvergencetest.py
              CMD_ARGS test_1pbox_network1d3d test_1pfv_network1d3d.input -Problem.Name test_1pbox_network1d3d)

# dim < dimWorld tests with Dune::Foamgrid<2,3>
dune_add_test(NAME test_1pbox_fracture2d3d
              SOURCES test_1pfv_fracture2d3d.cc
              COMPILE_DEFINITIONS TYPETAG=FractureBoxTypeTag
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/1pboxfracture2d3d-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/1pbox_fracture2d3d-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pbox_fracture2d3d test_1pfv_fracture2d3d.input -Problem.Name 1pbox_fracture2d3d")

dune_add_test(NAME test_1pcctpfa_fracture2d3d
              SOURCES test_1pfv_fracture2d3d.cc
              COMPILE_DEFINITIONS TYPETAG=FractureCCTpfaTypeTag
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/1pccfracture2d3d-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/1pcctpfa_fracture2d3d-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pcctpfa_fracture2d3d test_1pfv_fracture2d3d.input -Problem.Name 1pcctpfa_fracture2d3d")

dune_add_test(NAME test_1pccmpfa_fracture2d3d
              SOURCES test_1pfv_fracture2d3d.cc
              COMPILE_DEFINITIONS TYPETAG=FractureCCMpfaTypeTag
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/references/1pccfracture2d3d-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/1pccmpfa_fracture2d3d-00001.vtu
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pccmpfa_fracture2d3d test_1pfv_fracture2d3d.input -Problem.Name 1pccmpfa_fracture2d3d")

#install sources
install(FILES
1pniconductionproblem.hh
1pniconvectionproblem.hh
1pnispatialparams.hh
1ptestproblem.hh
1ptestspatialparams.hh
fractureproblem.hh
fracturespatialparams.hh
tubesproblem.hh
tubesspatialparams.hh
test_1pfv.cc
test_1pnifv.cc
test_1pfv_fracture2d3d.cc
test_1pfv_network1d3d.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/1p)
