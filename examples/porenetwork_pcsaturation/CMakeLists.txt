# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME example_pnm_pcsaturation
               LABELS porenetwork example
               SOURCES main.cc
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-foamgrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS --script fuzzyData --delimiter " "
                        --files ${CMAKE_SOURCE_DIR}/test/references/example_porenetwork_pcsaturation-reference.txt
                                ${CMAKE_CURRENT_BINARY_DIR}/pcsaturation_pnm_pc-s-curve.txt
                        --command "${CMAKE_CURRENT_BINARY_DIR}/example_pnm_pcsaturation params.input -Problem.PlotPcS false")
