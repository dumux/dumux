# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
add_input_file_links()

dumux_add_test(SOURCES test_decomposition.cc
               LABELS unit mortar
               CMAKE_GUARD "dune-foamgrid_FOUND")

dumux_add_test(SOURCES test_projectors.cc
               LABELS unit)

dumux_add_test(NAME test_trace_tpfa
               SOURCES test_trace.cc
               COMPILE_DEFINITIONS TYPETAG=TraceTpfa
               CMD_ARGS test_trace.input
               LABELS unit)
