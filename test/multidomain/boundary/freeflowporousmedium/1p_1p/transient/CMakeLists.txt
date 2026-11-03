# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_transient
              LABELS multidomain multidomain_boundary stokesdarcy 1p navierstokes
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_ff1p_pm1p_transient_freeflow-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_transient_freeflow-00019.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_ff1p_pm1p_transient_darcy-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_transient_darcy-00019.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_transient params.input
                                   -Vtk.OutputName test_md_boundary_ff1p_pm1p_transient")
