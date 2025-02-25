# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input solid.msh)
dumux_add_test(
    NAME test_saint_venant_kirchhoff_dynamic
    LABELS geomechanics elastic
    SOURCES main.cc
    CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
    COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
    CMD_ARGS  --script fuzzy
              --files ${CMAKE_SOURCE_DIR}/test/references/test_saint_venant_kirchhoff_dynamic-reference.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/test_saint_venant_kirchhoff_dynamic-00002.vtu
              --command "${CMAKE_CURRENT_BINARY_DIR}/test_saint_venant_kirchhoff_dynamic params.input -TimeLoop.TEnd 1 -VTKOutput.Every 50"
)
