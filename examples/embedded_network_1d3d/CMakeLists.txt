# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# link the parameter file and the grid file to the build directory
dune_symlink_to_source_files(FILES params.input network.dgf plot.py)

# add a new executable
# this could simply be
# dumux_add_test(
#     NAME example_embedded_network_1d3d
#     SOURCES main.cc
#)
# but we also do a regression test making sure that the example
# keep working when things in Dumux change
dumux_add_test(
    NAME example_embedded_network_1d3d
    SOURCES main.cc
    LABELS multidomain embedded example
    CMAKE_GUARD "( HAVE_UMFPACK AND dune-foamgrid_FOUND )"
    COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
    CMD_ARGS --script fuzzyData --delimiter " "
             --files ${CMAKE_SOURCE_DIR}/test/references/example_embedded_network_1d3d-reference.dat
                     ${CMAKE_CURRENT_BINARY_DIR}/clearance_tracer_amounts.dat
                     --command "${CMAKE_CURRENT_BINARY_DIR}/example_embedded_network_1d3d"
)
