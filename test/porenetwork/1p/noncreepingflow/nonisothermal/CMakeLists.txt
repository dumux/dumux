# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_pnm_1pni_noncreeping_flow
              SOURCES main.cc
              LABELS porenetwork
              COMPILE_DEFINITIONS ISOTHERMAL=0
              CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_1pni_noncreepingflow-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1pni_noncreepingflow-00014.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1pni_noncreeping_flow")
