dune_symlink_to_source_files(FILES "params.input" "cut.msh")

dumux_add_test(NAME test_richards_box_instationary
               LABELS porousmediumflow richards
               SOURCES main.cc
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_richards_box_instationary
               CMD_ARGS params.input -Problem.Name test_richards_box_instationary -SpatialParams.SoilType Clay)
