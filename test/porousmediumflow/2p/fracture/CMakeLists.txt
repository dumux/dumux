dune_symlink_to_source_files(FILES "params.input" "grids" "dt")

# tests without gravity
dumux_add_test(NAME test_2p_fracture_box
              LABELS porousmediumflow 2p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=FractureBox
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_fracture_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_box-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_box params.input -Problem.Name test_2p_fracture_box -TimeLoop.TimeStepSizeFile dt/test_2p_fracture_box_dt-reference.dat")

dumux_add_test(NAME test_2p_fracture_tpfa
              LABELS porousmediumflow 2p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=FractureCCTpfa
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_fracture_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_tpfa-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_tpfa params.input -Problem.Name test_2p_fracture_tpfa -TimeLoop.TimeStepSizeFile dt/test_2p_fracture_tpfa_dt-reference.dat")

dumux_add_test(NAME test_2p_fracture_mpfa
              LABELS porousmediumflow 2p
              SOURCES main.cc
              TIMEOUT 1500
              COMPILE_DEFINITIONS TYPETAG=FractureCCMpfa
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_fracture_mpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_mpfa-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_mpfa params.input -Problem.Name test_2p_fracture_mpfa -TimeLoop.TimeStepSizeFile dt/test_2p_fracture_mpfa_dt-reference.dat")

# tests with gravity
dumux_add_test(NAME test_2p_fracture_gravity_box
              TARGET test_2p_fracture_box
              LABELS porousmediumflow 2p
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_fracture_gravity_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_gravity_box-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_box params.input -Problem.Name test_2p_fracture_gravity_box -Problem.EnableGravity true -Newton.MaxRelativeShift 1e-12 -TimeLoop.TimeStepSizeFile dt/test_2p_fracture_gravity_box_dt-reference.dat")

dumux_add_test(NAME test_2p_fracture_gravity_tpfa
              TARGET test_2p_fracture_tpfa
              LABELS porousmediumflow 2p
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_fracture_gravity_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_gravity_tpfa-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_tpfa params.input -Problem.Name test_2p_fracture_gravity_tpfa -Problem.EnableGravity true -TimeLoop.TimeStepSizeFile dt/test_2p_fracture_gravity_tpfa_dt-reference.dat")

dumux_add_test(NAME test_2p_fracture_gravity_mpfa
              TARGET test_2p_fracture_mpfa
              LABELS porousmediumflow 2p
              TIMEOUT 1500
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_fracture_gravity_mpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_gravity_mpfa-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_fracture_mpfa params.input -Problem.Name test_2p_fracture_gravity_mpfa -Problem.EnableGravity true -TimeLoop.TimeStepSizeFile dt/test_2p_fracture_gravity_mpfa_dt-reference.dat")
