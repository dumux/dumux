# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input solid.msh)
dumux_add_test(
    NAME test_saint_venant_kirchhoff_stationary
    LABELS geomechanics elastic
    SOURCES main.cc
    CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
    COMMAND ${DUMUX_RUNTEST}
    CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_saint_venant_kirchhoff_stationary params.input -Problem.Name csm1"
)
