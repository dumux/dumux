# We execute Python tests in the binary directory so possible output files
# do not get dumped into the source tree

if(${dune-common_VERSION} VERSION_LESS 2.8)
  dune_python_add_test(NAME test_python_gridgeometry
                       COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_gridgeometry.py
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       LABELS python unit)

  dune_python_add_test(NAME test_python_fvproblem
                       COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_fvproblem.py
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       LABELS python unit)

  dune_python_add_test(NAME test_python_explicit_transport_cctpfa
                       LABELS python
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                               --script fuzzy
                               --files ${CMAKE_SOURCE_DIR}/test/references/test_python_explicit-transport-cctpfa-100.vtu
                                       ${CMAKE_CURRENT_BINARY_DIR}/finitevolume-solution-100.vtu
                               --command "${CMAKE_CURRENT_SOURCE_DIR}/test_explicit_transport_cctpfa.py")
else()
  dune_python_add_test(NAME test_python_gridgeometry
                       SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test_gridgeometry.py
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       LABELS python unit)

  dune_python_add_test(NAME test_python_fvproblem
                       SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test_fvproblem.py
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       LABELS python unit)

  dune_python_add_test(NAME test_python_explicit_transport_cctpfa
                       LABELS python
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       SCRIPT ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                              --script fuzzy
                              --files ${CMAKE_SOURCE_DIR}/test/references/test_python_explicit-transport-cctpfa-100.vtu
                                      ${CMAKE_CURRENT_BINARY_DIR}/finitevolume-solution-100.vtu
                              --command "${CMAKE_CURRENT_SOURCE_DIR}/test_explicit_transport_cctpfa.py")

  dune_python_add_test(NAME test_python_1p_incompressible_box_numdiff
                       LABELS python
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       SCRIPT ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                                --script fuzzy
                                --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_box-reference.vtu
                                        ${CMAKE_CURRENT_BINARY_DIR}/test_1p_box_numeric-00000.vtu
                                --command "${CMAKE_CURRENT_SOURCE_DIR}/test_1p.py box numeric")

  dune_python_add_test(NAME test_python_1p_incompressible_box_anadiff
                       LABELS python
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       SCRIPT ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                                --script fuzzy
                                --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_box-reference.vtu
                                        ${CMAKE_CURRENT_BINARY_DIR}/test_1p_box_analytic-00000.vtu
                                --command "${CMAKE_CURRENT_SOURCE_DIR}/test_1p.py box analytic")

  dune_python_add_test(NAME test_python_1p_incompressible_cctpfa_numdiff
                       LABELS python
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       SCRIPT ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                                --script fuzzy
                                --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_cc-reference.vtu
                                        ${CMAKE_CURRENT_BINARY_DIR}/test_1p_cctpfa_numeric-00000.vtu
                                --command "${CMAKE_CURRENT_SOURCE_DIR}/test_1p.py cctpfa numeric")

  dune_python_add_test(NAME test_python_1p_incompressible_cctpfa_anadiff
                       LABELS python
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       SCRIPT ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                                --script fuzzy
                                --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_cc-reference.vtu
                                        ${CMAKE_CURRENT_BINARY_DIR}/test_1p_cctpfa_analytic-00000.vtu
                                --command "${CMAKE_CURRENT_SOURCE_DIR}/test_1p.py cctpfa analytic")

endif()
