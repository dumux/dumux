# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# We execute Python tests in the binary directory so possible output files
# do not get dumped into the source tree

dune_python_add_test(NAME test_python_gridgeometry
                     SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test_gridgeometry.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                     LABELS python unit)
set_tests_properties(test_python_gridgeometry PROPERTIES TIMEOUT 1200)

dune_python_add_test(NAME test_python_fvproblem
                     SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test_fvproblem.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                     LABELS python unit)
set_tests_properties(test_python_fvproblem PROPERTIES TIMEOUT 1200)

dune_python_add_test(NAME test_python_explicit_transport_cctpfa
                     LABELS python
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                     SCRIPT ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                              --script fuzzy
                              --files ${CMAKE_SOURCE_DIR}/test/references/test_python_explicit-transport-cctpfa-100-reference.vtu
                                      ${CMAKE_CURRENT_BINARY_DIR}/finitevolume-solution-100.vtu
                              --command "${CMAKE_CURRENT_SOURCE_DIR}/test_explicit_transport_cctpfa.py")
set_tests_properties(test_python_explicit_transport_cctpfa PROPERTIES TIMEOUT 1200)

dune_python_add_test(NAME test_python_1p_incompressible_box_numdiff
                     LABELS python
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                     SCRIPT ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                                --script fuzzy
                                --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_box-reference.vtu
                                        ${CMAKE_CURRENT_BINARY_DIR}/test_1p_box_numeric-00000.vtu
                                --command "${CMAKE_CURRENT_SOURCE_DIR}/test_1p.py -DiscMethod box -DiffMethod numeric")
set_tests_properties(test_python_1p_incompressible_box_numdiff PROPERTIES TIMEOUT 1200)

dune_python_add_test(NAME test_python_1p_incompressible_box_anadiff
                     LABELS python
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                     SCRIPT ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                                --script fuzzy
                                --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_box-reference.vtu
                                        ${CMAKE_CURRENT_BINARY_DIR}/test_1p_box_analytic-00000.vtu
                                --command "${CMAKE_CURRENT_SOURCE_DIR}/test_1p.py -DiscMethod box -DiffMethod analytic")
set_tests_properties(test_python_1p_incompressible_box_anadiff PROPERTIES TIMEOUT 1200)

dune_python_add_test(NAME test_python_1p_incompressible_cctpfa_numdiff
                     LABELS python
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                     SCRIPT ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                                --script fuzzy
                                --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_cc-reference.vtu
                                        ${CMAKE_CURRENT_BINARY_DIR}/test_1p_cctpfa_numeric-00000.vtu
                                --command "${CMAKE_CURRENT_SOURCE_DIR}/test_1p.py -DiscMethod cctpfa -DiffMethod numeric")
set_tests_properties(test_python_1p_incompressible_cctpfa_numdiff PROPERTIES TIMEOUT 1200)

dune_python_add_test(NAME test_python_1p_incompressible_cctpfa_anadiff
                     LABELS python
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                     SCRIPT ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                                --script fuzzy
                                --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_cc-reference.vtu
                                        ${CMAKE_CURRENT_BINARY_DIR}/test_1p_cctpfa_analytic-00000.vtu
                                --command "${CMAKE_CURRENT_SOURCE_DIR}/test_1p.py -DiscMethod cctpfa -DiffMethod analytic")
set_tests_properties(test_python_1p_incompressible_cctpfa_anadiff PROPERTIES TIMEOUT 1200)
