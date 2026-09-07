dune_symlink_to_source_files(FILES "params_2d.input" "params_3d.input")

dumux_add_test(NAME test_2pve_columnmapping
               SOURCES test_columnmapping.cc
               LABELS unit porousmediumflow 2pve)

dumux_add_test(NAME test_2pve_quantityreconstruction
               SOURCES test_quantityreconstruction.cc
               LABELS unit porousmediumflow 2pve)

dumux_add_test(NAME test_2pve_immiscible_tpfa_2d
               SOURCES main.cc
               COMPILE_DEFINITIONS GRIDTYPE=Dune::YaspGrid<2,Dune::EquidistantOffsetCoordinates<double,2>>
               LABELS porousmediumflow 2pve
               TIMEOUT 3600
               COMMAND ${CMAKE_SOURCE_DIR}/../dumux/bin/testing/dumux_runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_2pve_coarse_2d_reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/2pve_lenses_2d_40x0x18-00067.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pve_immiscible_tpfa_2d params_2d.input -Problem.Name 2pve_lenses_2d_ -Grid.LowerLeft '0 0' -Grid.UpperRight '250 55' -Grid.Cells '40 18' -TimeLoop.TEnd 4e5 -TimeLoop.MaxTimeStepSize 10000 -BoundaryConditions.InjectionRate -5.835e-4")

dumux_add_test(NAME test_2pve_immiscible_tpfa_3d
               SOURCES main.cc
               COMPILE_DEFINITIONS GRIDTYPE=Dune::YaspGrid<3,Dune::EquidistantOffsetCoordinates<double,3>>
               LABELS porousmediumflow 2pve
               TIMEOUT 3600
               COMMAND ${CMAKE_SOURCE_DIR}/../dumux/bin/testing/dumux_runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_2pve_coarse_3d_reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/2pve_lenses_3d_15x8x3-00226.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pve_immiscible_tpfa_3d params_3d.input -Problem.Name 2pve_lenses_3d_ -Grid.LowerLeft '0 0 0' -Grid.UpperRight '250 100 25' -Grid.Cells '15 8 3' -TimeLoop.TEnd 2e6 -TimeLoop.MaxTimeStepSize 10000 -BoundaryConditions.InjectionRate -1.45875e-3")
