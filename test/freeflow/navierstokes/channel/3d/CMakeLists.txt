dune_add_test(NAME test_ff_stokes_channel_3d
              SOURCES main.cc
              COMPILE_DEFINITIONS DIM_3D=1
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel_3d-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d params.input
                             -Problem.Name test_ff_stokes_channel_3d"
                             --zeroThreshold {"velocity_liq \(m/s\)":1e-12})

dune_add_test(NAME test_ff_stokes_channel_pseudo3d
              SOURCES main.cc
              COMPILE_DEFINITIONS DIM_3D=0
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel_pseudo3d-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_pseudo3d-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_pseudo3d params_pseudo.input
                             -Problem.Name test_ff_stokes_channel_pseudo3d"
                             --zeroThreshold {"velocity_liq \(m/s\)":1e-12})


dune_symlink_to_source_files(FILES "params.input" "params_pseudo.input")
