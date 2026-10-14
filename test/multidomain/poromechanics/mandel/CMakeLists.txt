# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(analyticalsolution)

dune_symlink_to_source_files(FILES params.input plot_results.py convergencetest.py)

dumux_add_test(NAME test_md_mandel
              LABELS multidomain poromechanics 1p poroelastic
              SOURCES test_mandel.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_md_mandel_flow-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_md_mandel_flow-00002.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/test_md_mandel_mech-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_md_mandel_mech-00002.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_mandel params.input
                                  -TimeLoop.TEnd 50 -Vtk.OutputName test_md_mandel")

dumux_add_test(NAME test_md_mandel_compressible_grains
              LABELS multidomain poromechanics 1p poroelastic
              TARGET test_md_mandel
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_md_mandel_compressible_grains_flow-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_md_mandel_compressible_grains_flow-00002.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/test_md_mandel_compressible_grains_mech-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_md_mandel_compressible_grains_mech-00002.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_mandel params.input
                                  -TimeLoop.TEnd 50 -Problem.BiotCoefficient 0.7 -Problem.FluidBulkModulus 2.2e9
                                  -Vtk.OutputName test_md_mandel_compressible_grains")

dumux_add_test(NAME test_md_mandel_convergence
              LABELS multidomain poromechanics 1p poroelastic
              TARGET test_md_mandel
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/convergencetest.py
              CMD_ARGS test_md_mandel -Grid.Cells "20 2" -TimeLoop.Scheme DIRK3
                       -TimeLoop.TEnd 1000 -TimeLoop.CheckPoints 1000 -Vtk.OutputName test_md_mandel_convergence)
