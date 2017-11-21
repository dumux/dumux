dune_symlink_to_source_files(FILES "test_fracture.input" "grids")

dune_add_test(NAME test_2p_fracture_box
              SOURCES test_2p_fracture_fv.cc
              COMPILE_DEFINITIONS TYPETAG=FractureBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/fracturebox2p-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/fracturebox-00023.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_box test_fracture.input -Problem.Name fracturebox")

dune_add_test(NAME test_2p_fracture_tpfa
              SOURCES test_2p_fracture_fv.cc
              COMPILE_DEFINITIONS TYPETAG=FractureCCTpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/fracturecc2p-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/fracturetpfa-00029.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_tpfa test_fracture.input -Problem.Name fracturetpfa")

dune_add_test(NAME test_2p_fracture_mpfa
             SOURCES test_2p_fracture_fv.cc
             COMPILE_DEFINITIONS TYPETAG=FractureCCMpfaTypeTag
             COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
             CMD_ARGS --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/references/fractureccmpfa2p-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/fracturempfa-00031.vtu
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_mpfa test_fracture.input -Problem.Name fracturempfa")

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
problem.hh
spatialparams.hh
test_2p_fracture_box.cc
test_2p_fracture_cc.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2p/fracture)
