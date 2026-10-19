# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input hertz.geo run_benchmark.py)

dumux_add_test(NAME test_contact_hertz
               LABELS solidmechanics
               SOURCES main.cc
               CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK AND GMSH_FOUND )"
               COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/hertztest.py
               CMD_ARGS ${GMSH_EXECUTABLE} test_contact_hertz params.input -Problem.Approaches "0.0025 0.005")
