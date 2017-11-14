add_subdirectory(pointsources)
add_subdirectory(incompressible)
add_subdirectory(compressible)

add_input_file_links()
add_gstat_file_links()

dune_symlink_to_source_files(FILES "tubesconvergencetest.py")

# isothermal tests
dune_add_test(NAME test_1pcctpfa
              SOURCES test_1pfv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTestCCTpfaProblem
              COMMAND python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                      --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/1ptestcctpfa-00001.vtu
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pcctpfa test_1pfv.input -Problem.Name 1ptestcctpfa")

dune_add_test(NAME test_1pccmpfa
              SOURCES test_1pfv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTestCCMpfaProblem
              COMMAND python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                      --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/1ptestccmpfa-00001.vtu
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pccmpfa test_1pfv.input -Problem.Name 1ptestccmpfa")

dune_add_test(NAME test_1pbox
              SOURCES test_1pfv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTestBoxProblem
              COMMAND python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                      --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/references/1ptestbox-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/1ptestbox-00001.vtu
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pbox test_1pfv.input -Problem.Name 1ptestbox")

# add_dumux_test(test_cc1pwithgstat test_cc1pwithgstat test_cc1pwithgstat.cc
#                ./test_cc1pwithgstat)
#
# # non-isothermal tests
# add_dumux_test(test_box1pniconduction test_box1pniconduction test_box1pniconduction.cc
#                python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
#                  --script fuzzy
#                  --files ${CMAKE_SOURCE_DIR}/test/references/1pniboxconduction-reference.vtu
#                          ${CMAKE_CURRENT_BINARY_DIR}/box1pniconduction-00006.vtu
#                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1pniconduction"
#                  --zeroThreshold {"velocity_H2O":1e-8})
#
# add_dumux_test(test_box1pniconvection test_box1pniconvection test_box1pniconvection.cc
#                python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
#                  --script fuzzy
#                  --files ${CMAKE_SOURCE_DIR}/test/references/1pniboxconvection-reference.vtu
#                          ${CMAKE_CURRENT_BINARY_DIR}/box1pniconvection-00010.vtu
#                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1pniconvection"
#                  --zeroThreshold {"velocity":1e-15})
#
# add_dumux_test(test_cc1pniconduction test_cc1pniconduction test_cc1pniconduction.cc
#                python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
#                  --script fuzzy
#                  --files ${CMAKE_SOURCE_DIR}/test/references/1pniccconduction-reference.vtu
#                          ${CMAKE_CURRENT_BINARY_DIR}/cc1pniconduction-00006.vtu
#                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1pniconduction"
#                  --zeroThreshold {"velocity":1e-8})
#
# add_dumux_test(test_cc1pniconvection test_cc1pniconvection test_cc1pniconvection.cc
#                python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
#                  --script fuzzy
#                  --files ${CMAKE_SOURCE_DIR}/test/references/1pniccconvection-reference.vtu
#                          ${CMAKE_CURRENT_BINARY_DIR}/cc1pniconvection-00010.vtu
#                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1pniconvection")
#
# add_dumux_test(test_ccmpfa1pniconduction test_ccmpfa1pniconduction test_ccmpfa1pniconduction.cc
#                python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
#                  --script fuzzy
#                  --files ${CMAKE_SOURCE_DIR}/test/references/ccmpfa1pniconduction-reference.vtu
#                          ${CMAKE_CURRENT_BINARY_DIR}/ccmpfa1pniconduction-00006.vtu
#                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccmpfa1pniconduction")
#
# add_dumux_test(test_ccmpfa1pniconvection test_ccmpfa1pniconvection test_ccmpfa1pniconvection.cc
#                python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
#                  --script fuzzy
#                  --files ${CMAKE_SOURCE_DIR}/test/references/ccmpfa1pniconvection-reference.vtu
#                          ${CMAKE_CURRENT_BINARY_DIR}/ccmpfa1pniconvection-00010.vtu
#                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccmpfa1pniconvection")
#
# dim < dimWorld tests with foamgrid
dune_add_test(NAME test_1pcctpfa_network1d3d
              SOURCES test_1pfv_network1d3d.cc
              COMPILE_DEFINITIONS TYPETAG=TubesTestCCTpfaProblem
              COMMAND ./tubesconvergencetest.py test_1pcctpfa_network1d3d test_1pfv_network1d3d.input -Problem.Name test_1pcctpfa_network1d3d)

dune_add_test(NAME test_1pbox_network1d3d
              SOURCES test_1pfv_network1d3d.cc
              COMPILE_DEFINITIONS TYPETAG=TubesTestBoxProblem
              COMMAND ./tubesconvergencetest.py test_1pbox_network1d3d test_1pfv_network1d3d.input -Problem.Name test_1pbox_network1d3d)

# dune_add_test(SOURCES test_box1pnetwork1d3d.cc
#               NAME test_box1pnetwork1d3d
#               COMMAND ./tubesconvergencetest.py test_box1pnetwork1d3d)
#
# add_dumux_test(test_box1pfracture2d3d test_box1pfracture2d3d test_box1pfracture2d3d.cc
#                python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
#                  --script fuzzy
#                  --files ${CMAKE_SOURCE_DIR}/test/references/1pboxfracture2d3d-reference.vtu
#                          ${CMAKE_CURRENT_BINARY_DIR}/1pboxfracture2d3d-00001.vtu
#                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1pfracture2d3d")
#
# dune_add_test(SOURCES test_cc1pnetwork1d3d.cc
#               NAME test_cc1pnetwork1d3d
#               COMMAND ./tubesconvergencetest.py test_cc1pnetwork1d3d)
#
# add_dumux_test(test_cc1pfracture2d3d test_cc1pfracture2d3d test_cc1pfracture2d3d.cc
#                python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
#                  --script fuzzy
#                  --files ${CMAKE_SOURCE_DIR}/test/references/1pccfracture2d3d-reference.vtu
#                          ${CMAKE_CURRENT_BINARY_DIR}/1pccfracture2d3d-00001.vtu
#                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1pfracture2d3d")

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
test_box1pfracture2d3d.cc
test_box1pnetwork1d3d.cc
test_box1p.cc
test_box1pniconduction.cc
test_box1pniconvection.cc
test_box1pwithamg.cc
test_cc1pfracture2d3d.cc
test_cc1pnetwork1d3d.cc
test_cc1p.cc
test_cc1pniconduction.cc
test_cc1pniconvection.cc
test_cc1pwithamg.cc
test_ccmpfa1pniconduction.cc
test_ccmpfa1pniconvection.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/1p)
