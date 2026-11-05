# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input" "plot_convergence.py")

# injection rate of 120 m^3/day of Nordbotten and Celia (2006), gravity number 0.14
dumux_add_test(NAME test_2pve_radialinjection_tpfa
               LABELS porousmediumflow 2pve
               SOURCES main.cc
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_2pve_radialinjection_tpfa
               CMD_ARGS params.input)

# the same injected volume at a hundred times the injection rate, gravity number 0.0014
dumux_add_test(NAME test_2pve_radialinjection_smallgravity_tpfa
               TARGET test_2pve_radialinjection_tpfa
               LABELS porousmediumflow 2pve
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_2pve_radialinjection_tpfa
               CMD_ARGS params.input
                        -Problem.Name test_2pve_radialinjection_smallgravity_tpfa
                        -Problem.InjectionRate 1.388888889e-1
                        -TimeLoop.TEnd 8.64e6
                        -TimeLoop.MaxTimeStepSize 8.64e3
                        -Benchmark.MaxRelativeGasPlumeDistanceError 0.005)

# negligible gravity number (1.4e-4) and capillary fringe (entry pressure 1 Pa), maximum time step size refined with the radial cell size
dumux_add_test(NAME test_2pve_radialinjection_convergence
               TARGET test_2pve_radialinjection_tpfa
               LABELS porousmediumflow 2pve
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pve_radialinjection_tpfa params.input -Problem.Name test_2pve_radialinjection_convergence -Problem.InjectionRate 1.388888889 -TimeLoop.TEnd 8.64e5 -TimeLoop.DtInitial 0.01 -SpatialParams.BrooksCoreyPcEntry 1 -Benchmark.MaxRelativeGasPlumeDistanceError 1"
                        --levels "400 864" "800 432"
                        --level-args "-Grid.Cells '{0} 30' -TimeLoop.MaxTimeStepSize {1}"
                        --errors test_2pve_radialinjection_convergence-errors.csv
                        --check "L1Rel(zp_eff) mean 0.8 1.2")
