dune_symlink_to_source_files(FILES "params.input" "periodic.dgf")

# test using periodic boundary conditions
dumux_add_test(NAME test_1p_periodic_tpfa
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleTpfa FVGEOMCACHING=false
              CMAKE_GUARD dune-spgrid_FOUND
              COMPILE_ONLY)

dumux_add_test(NAME test_1p_periodic_tpfa_sequential
              TARGET test_1p_periodic_tpfa
              LABELS porousmediumflow 1p
              CMAKE_GUARD dune-spgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_tpfa_periodic-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_periodic_tpfa_sequential-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_periodic_tpfa params.input -Problem.Name test_1p_periodic_tpfa_sequential")

dumux_add_test(NAME test_1p_periodic_tpfa_parallel
              TARGET test_1p_periodic_tpfa
              LABELS porousmediumflow 1p parallel
              CMAKE_GUARD dune-spgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --zeroThreshold {"process rank":100}
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_tpfa_periodic-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_1p_periodic_tpfa_parallel-00001.pvtu
                        --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_1p_periodic_tpfa params.input -Problem.Name test_1p_periodic_tpfa_parallel")

dumux_add_test(NAME test_1p_periodic_tpfa_caching
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleTpfa FVGEOMCACHING=true
              CMAKE_GUARD dune-spgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_tpfa_periodic-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_periodic_tpfa_caching-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_periodic_tpfa params.input -Problem.Name test_1p_periodic_tpfa_caching")

dumux_add_test(NAME test_1p_periodic_box
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleBox
              CMAKE_GUARD dune-spgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_box_periodic-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_periodic_box-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_periodic_box params.input -Problem.Name test_1p_periodic_box")
