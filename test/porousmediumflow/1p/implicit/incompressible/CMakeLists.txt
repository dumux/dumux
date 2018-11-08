dune_symlink_to_source_files(FILES "params.input")

# using tpfa and analytical Jacobian
dune_add_test(NAME test_1p_incompressible_tpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleTpfa NUMDIFFMETHOD=DiffMethod::analytic
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_cc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_tpfa-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_tpfa params.input -Problem.Name test_1p_incompressible_tpfa")

# using tpfa, analytical Jacobian and quad precision
dune_add_test(NAME test_1p_incompressible_tpfa_quad
              SOURCES main.cc
              CMAKE_GUARD HAVE_QUAD
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleTpfaQuad NUMDIFFMETHOD=DiffMethod::analytic
              LINK_LIBRARIES ${QUADMATH_LIBRARIES}
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_cc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_tpfa_quad-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_tpfa_quad params.input -Problem.Name test_1p_incompressible_tpfa_quad")

# using mpfa and analytical Jacobian
dune_add_test(NAME test_1p_incompressible_mpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleMpfa NUMDIFFMETHOD=DiffMethod::analytic
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_cc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_mpfa-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_mpfa params.input -Problem.Name test_1p_incompressible_mpfa")

# using box and analytical Jacobian
dune_add_test(NAME test_1p_incompressible_box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleBox NUMDIFFMETHOD=DiffMethod::analytic
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_box-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_box params.input -Problem.Name test_1p_incompressible_box")


# using tpfa and numeric differentiation
dune_add_test(NAME test_1p_incompressible_tpfa_numdiff
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleTpfa NUMDIFFMETHOD=DiffMethod::numeric
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_cc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_tpfa_numdiff-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_tpfa_numdiff params.input
                                                                                  -Problem.Name test_1p_incompressible_tpfa_numdiff
                                                                                  -Assembly.NumericDifference.PriVarMagnitude 1e5")

# using box and numeric differentiation
dune_add_test(NAME test_1p_incompressible_box_numdiff
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleBox NUMDIFFMETHOD=DiffMethod::numeric
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_box_numdiff-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_box_numdiff params.input
                                                                                  -Problem.Name test_1p_incompressible_box_numdiff
                                                                                  -Assembly.NumericDifference.PriVarMagnitude 1e5")

set(CMAKE_BUILD_TYPE Release)
