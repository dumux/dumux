# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()
dune_symlink_to_source_files(FILES voronoi.pbm)

dumux_add_test(NAME test_periodic_porescale
               SOURCES main.cc
               LABELS porescale periodic freeflow
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-spgrid_FOUND AND dune-subgrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_periodic_subgrid-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_periodic_subgrid-00000.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_ff_periodic_subgrid_tiled-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_periodic_subgrid_tiled-00000.vtu
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_periodic_porescale params.input")
