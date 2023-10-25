# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_ff_stokes2c_densitydrivenflow
              LABELS freeflow navierstokesnc navierstokes
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes2c_densitydrivenflow-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c_densitydrivenflow-00021.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c_densitydrivenflow params.input
                             -Problem.Name test_ff_stokes2c_densitydrivenflow")
