# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(chrono)
add_subdirectory(container)
add_subdirectory(format)
add_subdirectory(gnuplotinterface)
add_subdirectory(gridmanager)
add_subdirectory(inputdata)
add_subdirectory(rasterimagereader)
add_subdirectory(vtk)

dumux_add_test(
    NAME test_grid_writer
    SOURCES test_grid_writer.cc
    CMAKE_GUARD DUMUX_HAVE_GRIDFORMAT
)

dumux_add_test(
    NAME test_grid_writer_parallel
    SOURCES test_grid_writer_parallel.cc
    CMAKE_GUARD DUMUX_HAVE_GRIDFORMAT
    MPI_RANKS 1 2
    TIMEOUT 100
)
