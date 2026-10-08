# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input convergencetest.py)

dumux_add_test(NAME test_ff_stokes_axisymmetric_pq1bubble
               LABELS freeflow navierstokes
               SOURCES main.cc
               CMAKE_GUARD "( dune-uggrid_FOUND AND HAVE_UMFPACK )"
               COMPILE_DEFINITIONS TYPETAG=AxisymmetricStokesPQ1Bubble
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/convergencetest.py
               CMD_ARGS ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_axisymmetric_pq1bubble 1.8 -Grid.CellType Simplex)

dumux_add_test(NAME test_ff_stokes_axisymmetric_pq1bubble_integral
               LABELS freeflow navierstokes
               SOURCES main.cc
               CMAKE_GUARD "( dune-uggrid_FOUND AND HAVE_UMFPACK )"
               COMPILE_DEFINITIONS TYPETAG=AxisymmetricStokesPQ1Bubble
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/convergencetest.py
               CMD_ARGS ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_axisymmetric_pq1bubble_integral 1.8 -Grid.CellType Simplex)

dumux_add_test(NAME test_ff_stokes_axisymmetric_pq1bubblehybrid
               LABELS freeflow navierstokes
               SOURCES main.cc
               CMAKE_GUARD "( dune-uggrid_FOUND AND HAVE_UMFPACK )"
               COMPILE_DEFINITIONS TYPETAG=AxisymmetricStokesPQ1BubbleHybrid
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/convergencetest.py
               CMD_ARGS ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_axisymmetric_pq1bubblehybrid 1.8 -Grid.CellType Cube)

dumux_add_test(NAME test_ff_stokes_axisymmetric_pq2hybrid
               LABELS freeflow navierstokes
               SOURCES main.cc
               CMAKE_GUARD "( dune-uggrid_FOUND AND HAVE_UMFPACK )"
               COMPILE_DEFINITIONS TYPETAG=AxisymmetricStokesPQ2Hybrid
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/convergencetest.py
               CMD_ARGS ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_axisymmetric_pq2hybrid 1.8 -Grid.CellType Cube)
