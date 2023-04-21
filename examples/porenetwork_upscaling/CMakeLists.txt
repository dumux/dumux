# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME example_pnm1p_upscaling
               LABELS porenetwork example
               SOURCES main.cc
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-foamgrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS       --script fuzzyData --delimiter " "
                              --files ${CMAKE_SOURCE_DIR}/test/references/example_porenetwork_upscaling_noncreeping_flow_X-dir_PermeabilityratioVsForchheimerNumber.dat
                                      ${CMAKE_CURRENT_BINARY_DIR}/X-dir-PermeabilityratioVsForchheimerNumber.dat
                                      ${CMAKE_SOURCE_DIR}/test/references/example_porenetwork_upscaling_noncreeping_flow_X-dir_InversePrmeabilityVsInertiaToViscousRatio.dat
                                      ${CMAKE_CURRENT_BINARY_DIR}/X-dir-InversePrmeabilityVsInertiaToViscousRatio.dat
                              --command "${CMAKE_CURRENT_BINARY_DIR}/example_pnm1p_upscaling params.input -Problem.AssumeCreepingFlow false -Problem.Directions 0 -Output.PlotConductivity false")

dumux_add_test(NAME example_pnm1p_creeping_flow_upscaling
               TARGET example_pnm1p_upscaling
               LABELS porenetwork example
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-foamgrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/example_pnm1p_upscaling  params.input -Problem.AssumeCreepingFlow true")
