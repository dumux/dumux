# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_2ptracer_lens_tpfa
              LABELS porousmediumflow tracer 2p
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2ptracer_lens_tpfa_transport-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2ptracer_lens_tpfa-00050.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/test_2ptracer_lens_tpfa_2p-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2ptracer_lens_tpfa_2p-00050.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2ptracer_lens_tpfa params.input -Problem.Name test_2ptracer_lens_tpfa")
