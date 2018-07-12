dune_symlink_to_source_files(FILES "test_1p.input" "periodic.dgf")

set(CMAKE_BUILD_TYPE Release)

# test using periodic boundary conditions
dune_add_test(NAME test_1p_periodic_tpfa
              SOURCES test_1pfv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleTpfa FVGEOMCACHING=false
              CMAKE_GUARD dune-spgrid_FOUND
              COMPILE_ONLY)

dune_add_test(NAME test_1p_periodic_tpfa_sequential
              TARGET test_1p_periodic_tpfa
              CMAKE_GUARD dune-spgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test-1p-tpfa-periodic-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_tpfa-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_periodic_tpfa test_1p.input -Problem.Name test_1p_tpfa")

dune_add_test(NAME test_1p_periodic_tpfa_parallel
              TARGET test_1p_periodic_tpfa
              CMAKE_GUARD dune-spgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --zeroThreshold {"process rank":100}
                        --files ${CMAKE_SOURCE_DIR}/test/references/test-1p-tpfa-periodic-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_1p_tpfa_parallel-00001.pvtu
                        --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_1p_periodic_tpfa test_1p.input -Problem.Name test_1p_tpfa_parallel")

dune_add_test(NAME test_1p_periodic_tpfa_caching
              SOURCES test_1pfv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleTpfa FVGEOMCACHING=true
              CMAKE_GUARD dune-spgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test-1p-tpfa-periodic-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_tpfa_caching-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_periodic_tpfa test_1p.input -Problem.Name test_1p_tpfa_caching")

#dune_add_test(NAME test_1p_periodic_box
#              SOURCES test_1pfv.cc
#              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleBox
#              CMAKE_GUARD dune-spgrid_FOUND
#              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
#              CMD_ARGS  --script fuzzy
#                        --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
#                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_box-00001.vtu
#                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_periodic_box test_1p.input -Problem.Name test_1p_box")
