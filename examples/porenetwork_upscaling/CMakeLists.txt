dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME example_pnm1p_upscaling
               LABELS porenetwork example
               SOURCES main.cc
               CMAKE_GUARD HAVE_UMFPACK dune-foamgrid_FOUND
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS       --script fuzzyData --delimiter " "
                              --files ${CMAKE_SOURCE_DIR}/test/references/example_porenetwork_upscaling_noncreeping_flow_X_direction-reference.dat
                                      ${CMAKE_CURRENT_BINARY_DIR}/X-dir.dat
                                      ${CMAKE_SOURCE_DIR}/test/references/example_porenetwork_upscaling_noncreeping_flow_Y_direction-reference.dat
                                      ${CMAKE_CURRENT_BINARY_DIR}/Y-dir.dat
                                      ${CMAKE_SOURCE_DIR}/test/references/example_porenetwork_upscaling_noncreeping_flow_Z_direction-reference.dat
                                      ${CMAKE_CURRENT_BINARY_DIR}/Z-dir.dat
                              --command "${CMAKE_CURRENT_BINARY_DIR}/example_pnm1p_upscaling params.input -Problem.AssumeCreepingFlow false")

dumux_add_test(NAME example_pnm1p_creeping_flow_upscaling
               TARGET example_pnm1p_upscaling
               CMAKE_GUARD HAVE_UMFPACK dune-foamgrid_FOUND
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/example_pnm1p_upscaling  params.input -Problem.AssumeCreepingFlow true")