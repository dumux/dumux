dune_symlink_to_source_files(FILES "params.input" "grids")

# tests without gravity
dune_add_test(NAME test_2p_fracture_box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=FractureBox
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_fracture_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_box-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_box params.input -Problem.Name test_2p_fracture_box")

dune_add_test(NAME test_2p_fracture_tpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=FractureCCTpfa
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_fracture_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_tpfa-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_tpfa params.input -Problem.Name test_2p_fracture_tpfa")

dune_add_test(NAME test_2p_fracture_mpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=FractureCCMpfa
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_fracture_mpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_mpfa-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_mpfa params.input -Problem.Name test_2p_fracture_mpfa")

# tests with gravity
dune_add_test(NAME test_2p_fracture_gravity_box
              TARGET test_2p_fracture_box
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_fracture_gravity_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_gravity_box-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_box params.input -Problem.Name test_2p_fracture_gravity_box -Problem.EnableGravity true")

dune_add_test(NAME test_2p_fracture_gravity_tpfa
              TARGET test_2p_fracture_tpfa
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_fracture_gravity_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_gravity_tpfa-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_tpfa params.input -Problem.Name test_2p_fracture_gravity_tpfa -Problem.EnableGravity true")

dune_add_test(NAME test_2p_fracture_gravity_mpfa
              TARGET test_2p_fracture_mpfa
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_fracture_gravity_mpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_gravity_mpfa-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_mpfa params.input -Problem.Name test_2p_fracture_gravity_mpfa -Problem.EnableGravity true")

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
problem.hh
spatialparams.hh
main.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2p/fracture)
