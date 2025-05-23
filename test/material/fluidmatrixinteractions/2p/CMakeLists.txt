# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(SOURCES test_thermalconductivity.cc
              LABELS unit material
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzyData --delimiter " "
                      --files ${CMAKE_SOURCE_DIR}/test/references/test_thermalconductivitysomerton-reference.dat
                              ${CMAKE_CURRENT_BINARY_DIR}/somerton_lambda_eff.dat
                              ${CMAKE_SOURCE_DIR}/test/references/test_thermalconductivityjohansen-reference.dat
                              ${CMAKE_CURRENT_BINARY_DIR}/johansen_lambda_eff.dat
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_thermalconductivity")

dumux_add_test(SOURCES test_material_2p_vangenuchten.cc
               LABELS unit material
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS --script fuzzyData --delimiter " "
                        --files  ${CMAKE_SOURCE_DIR}/test/references/test_pcsw_vangenuchten-reference.dat
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_pcsw_vangenuchten.dat
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_material_2p_vangenuchten")

dumux_add_test(SOURCES test_material_2p_brookscorey.cc
               LABELS unit material
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS --script fuzzyData --delimiter " "
                        --files  ${CMAKE_SOURCE_DIR}/test/references/test_pcsw_brookscorey-reference.dat
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_pcsw_brookscorey.dat
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_material_2p_brookscorey")


dumux_add_test(SOURCES test_material_2p_spline.cc
               LABELS unit material)

dumux_add_test(SOURCES test_material_2p_dataspline.cc
               LABELS unit material)

dune_symlink_to_source_files(FILES test_material_2p_spline.input test_material_2p_dataspline.input)
