# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input grid.msh)

dumux_add_test(
    NAME test_hyperelastic_thick_cylindrical_shell
    LABELS geomechanics elastic
    SOURCES main.cc
    CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_hyperelastic_thick_cylindrical_shell
    CMD_ARGS params.input
)
