if(${dune-common_VERSION} VERSION_LESS 2.8)
  dune_python_add_test(NAME test_python_gridgeometry
                       COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_py_gridgeometry.py
                       LABELS python unit)

  dune_python_add_test(NAME test_python_fvproblem
                       COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_py_fvproblem.py
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
  dune_python_add_test(NAME test_py_gridgeometry
                       SCRIPT test_py_gridgeometry.py
                       LABELS python unit)

  dune_python_add_test(NAME test_py_fvproblem
                       SCRIPT test_py_fvproblem.py
                       LABELS python unit)

  dune_python_add_test(NAME test_python_explicit_transport_cctpfa
                       LABELS python
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       SCRIPT ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                              --script fuzzy
                              --files ${CMAKE_SOURCE_DIR}/test/references/test_python_explicit-transport-cctpfa-100.vtu
                                      ${CMAKE_CURRENT_BINARY_DIR}/finitevolume-solution-100.vtu
                              --command "${CMAKE_CURRENT_SOURCE_DIR}/test_explicit_transport_cctpfa.py")

endif()
