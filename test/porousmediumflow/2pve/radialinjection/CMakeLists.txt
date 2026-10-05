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
                        -Benchmark.MaxRelativeInterfaceError 0.005)
