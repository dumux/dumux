dune_symlink_to_source_files(FILES "test_1p.input")

# using tpfa and analytical Jacobian
dune_add_test(NAME test_1p_incompressible_tpfa
              SOURCES test_1pfv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleTpfa NUMDIFFMETHOD=DiffMethod::analytic
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1ptestcctpfa-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_tpfa test_1p.input -Problem.Name 1ptestcctpfa")

# using tpfa, analytical Jacobian and quad precision
dune_add_test(NAME test_1p_incompressible_tpfa_quad
              SOURCES test_1pfv.cc
              CMAKE_GUARD HAVE_QUAD
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleTpfaQuad NUMDIFFMETHOD=DiffMethod::analytic
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1ptestcctpfaquad-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_tpfa_quad test_1p.input -Problem.Name 1ptestcctpfaquad")
target_link_libraries(test_1p_incompressible_tpfa_quad ${QUADMATH_LIBRARIES})

# using mpfa and analytical Jacobian
dune_add_test(NAME test_1p_incompressible_mpfa
              SOURCES test_1pfv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleMpfa NUMDIFFMETHOD=DiffMethod::analytic
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1ptestccmpfa-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_mpfa test_1p.input -Problem.Name 1ptestccmpfa")

# using box and analytical Jacobian
dune_add_test(NAME test_1p_incompressible_box
              SOURCES test_1pfv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleBox NUMDIFFMETHOD=DiffMethod::analytic
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1ptestbox-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1ptestbox-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_box test_1p.input -Problem.Name 1ptestbox")


# using tpfa and numeric differentiation
dune_add_test(NAME test_1p_incompressible_tpfa_numdiff
              SOURCES test_1pfv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleTpfa NUMDIFFMETHOD=DiffMethod::numeric
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1ptestcctpfa_numdiff-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_tpfa_numdiff test_1p.input
                                                                                  -Problem.Name 1ptestcctpfa_numdiff
                                                                                  -Assembly.NumericDifference.PriVarMagnitude 1e5")

# using box and numeric differentiation
dune_add_test(NAME test_1p_incompressible_box_numdiff
              SOURCES test_1pfv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleBox NUMDIFFMETHOD=DiffMethod::numeric
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1ptestbox-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1ptestbox_numdiff-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_box_numdiff test_1p.input
                                                                                  -Problem.Name 1ptestbox_numdiff
                                                                                  -Assembly.NumericDifference.PriVarMagnitude 1e5")

set(CMAKE_BUILD_TYPE Release)
