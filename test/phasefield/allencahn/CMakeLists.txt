# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_phasefield_allencahn_pq2
               LABELS phasefield allencahn
               SOURCES main.cc
               CMAKE_GUARD DUMUX_HAVE_GRIDFORMAT
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_phasefield_allencahn_pq2-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_phasefield_allencahn_pq2-00023.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_phasefield_allencahn_pq2 params.input
                        -Grid.Cells \" 10 10\" -Output.Every 1000")

# Quantitative curvature-driven-shrinking benchmark: parameters chosen so the
# equilibrium interface width (sqrt(2*SurfaceTension/EnergyScale) = 0.04) is
# well resolved by the grid (h=0.01), unlike the qualitative demo test above
# (whose narrower equilibrium width is under-resolved and shows numerically
# pinned, not physically accurate, interface motion -- see the workpackage
# notes). With this resolved width, the discrete solution's shrinking rate
# matches the sharp-interface prediction dR/dt=-Mobility*SurfaceTension/R to
# ~2% (checked in main.cc's "Shrinking-circle rate check").
dumux_add_test(NAME test_phasefield_allencahn_shrinkingcircle
               LABELS phasefield allencahn
               TARGET test_phasefield_allencahn_pq2
               CMAKE_GUARD DUMUX_HAVE_GRIDFORMAT
               CMD_ARGS params.input
                        -Problem.Name test_phasefield_allencahn_shrinkingcircle
                        -Problem.InitialRadius 0.2 -Problem.SurfaceTension 0.0004 -Problem.EnergyScale 0.5
                        -Problem.Mobility 50 -TimeLoop.TEnd 0.5 -TimeLoop.InitialTimeStepSize 1e-3s
                        -TimeLoop.MaxTimeStepSize 5e-3s -Output.Every 100000)
