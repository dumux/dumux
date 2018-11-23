dune_symlink_to_source_files(FILES "params.input" "grids")

# using box
add_executable(test_2p_fracture_box EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_2p_fracture_box PUBLIC TYPETAG=FractureBox)

# using tpfa
add_executable(test_2p_fracture_tpfa EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_2p_fracture_tpfa PUBLIC TYPETAG=FractureCCTpfa)

# using tpfa
add_executable(test_2p_fracture_mpfa EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_2p_fracture_mpfa PUBLIC TYPETAG=FractureCCMpfa)

# box test without gravity
dune_add_test(NAME test_2p_fracture_nogravity_box
              TARGET test_2p_fracture_box
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_fracture_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_box-00023.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_box params.input -Problem.Name test_2p_fracture_box")

dune_add_test(NAME test_2p_fracture_nogravity_tpfa
              TARGET test_2p_fracture_tpfa
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_fracture_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_tpfa-00029.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_tpfa params.input -Problem.Name test_2p_fracture_tpfa")

dune_add_test(NAME test_2p_fracture_nogravity_mpfa
              TARGET test_2p_fracture_mpfa
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_fracture_mpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_mpfa-00030.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_mpfa params.input -Problem.Name test_2p_fracture_mpfa")

# tests with gravity
dune_add_test(NAME test_2p_fracture_gravity_box
              TARGET test_2p_fracture_box
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_fracture_gravity_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_gravity_box-00028.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_box params.input -Problem.Name test_2p_fracture_gravity_box -Problem.EnableGravity true")

dune_add_test(NAME test_2p_fracture_gravity_tpfa
              TARGET test_2p_fracture_tpfa
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_fracture_gravity_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_gravity_tpfa-00042.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_tpfa params.input -Problem.Name test_2p_fracture_gravity_tpfa -Problem.EnableGravity true")

dune_add_test(NAME test_2p_fracture_gravity_mpfa
              TARGET test_2p_fracture_mpfa
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_fracture_gravity_mpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_gravity_mpfa-00042.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_mpfa params.input -Problem.Name test_2p_fracture_gravity_mpfa -Problem.EnableGravity true")

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
problem.hh
spatialparams.hh
main.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2p/fracture)
