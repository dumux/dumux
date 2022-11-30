# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(container)
add_subdirectory(format)
add_subdirectory(gnuplotinterface)
add_subdirectory(gridmanager)
add_subdirectory(inputdata)
add_subdirectory(rasterimagereader)
add_subdirectory(vtk)

dumux_add_test(
    NAME test_timeseries_writer
    SOURCES test_timeseries_writer.cc
    CMAKE_GUARD HAVE_GRIDFORMAT
)
