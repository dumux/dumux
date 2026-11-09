# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# additional macros
include(AddInputFileLinks)
include(DumuxDoxygen)
include(DumuxTestMacros)

find_package(Gnuplot QUIET)
set(DUMUX_HAVE_GNUPLOT ${GNUPLOT_FOUND})
find_package(Gstat QUIET)
find_package(Gmsh QUIET)
find_package(PTScotch QUIET)
# CGAL's configuration resets the policies and looks for Boost with the module that CMP0167 removes
set(CMAKE_POLICY_DEFAULT_CMP0167 NEW)
find_package(CGAL QUIET)
unset(CMAKE_POLICY_DEFAULT_CMP0167)
include(AddPTScotchFlags)
