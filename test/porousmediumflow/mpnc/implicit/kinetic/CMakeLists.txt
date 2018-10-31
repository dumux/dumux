add_input_file_links()

# build target for the full kinetic test problem
dune_add_test(NAME test_mpnc_kinetic_box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=EvaporationAtmosphereBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_mpnc_kinetic_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_kinetic_box-00011.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_kinetic_box params.input -Problem.Name test_mpnc_kinetic_box")

#install sources
install(FILES
problem.hh
spatialparams.hh
plotoverline2d.hh
main.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/mpnc/kinetic)
