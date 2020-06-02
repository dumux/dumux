add_input_file_links()
dune_symlink_to_source_files(FILES grids)

# build target for the CO2 test problem
# Ignore the porosity for all box models since it is defined element-wise in these test
# but the default 2p2c implementation outputs porosity per vertex.
# Depending on the order of the elements, the porosity would differ in these cases.
add_executable(test_co2_box EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_co2_box PUBLIC TYPETAG=HeterogeneousBox )

add_executable(test_co2_tpfa EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_co2_tpfa PUBLIC TYPETAG=HeterogeneousCCTpfa )

# sequential box
dumux_add_test(NAME test_co2_box_sequential
               TARGET test_co2_box
               LABELS porousmediumflow co2
               CMAKE_GUARD "( dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_co2_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_co2_box-00028.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_co2_box params.input -Problem.Name test_co2_box"
                        --zeroThreshold {"porosity":1})

# sequential tpfa
dumux_add_test(NAME test_co2_tpfa_sequential
               LABELS porousmediumflow co2
               TARGET test_co2_tpfa
               COMPILE_DEFINITIONS TYPETAG=HeterogeneousCCTpfa
               CMAKE_GUARD "( dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_co2_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_co2_tpfa-00028.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_co2_tpfa params.input -Problem.Name test_co2_tpfa")

# parallel box
dumux_add_test(NAME test_co2_box_parallel
               TARGET test_co2_box
               LABELS porousmediumflow co2 parallel
               TIMEOUT 1200
               CMAKE_GUARD "( dune-alugrid_FOUND AND MPI_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy --zeroThreshold {"process rank":100,"porosity":1}
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_co2_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_co2_box_parallel-00028.pvtu
                        --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_co2_box params.input -Problem.Name test_co2_box_parallel")

# parallel tpfa
dumux_add_test(NAME test_co2_tpfa_parallel
               TARGET test_co2_tpfa
               LABELS porousmediumflow co2 parallel
               TIMEOUT 1200
               CMAKE_GUARD "( dune-alugrid_FOUND AND MPI_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy --zeroThreshold {"process rank":100}
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_co2_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_co2_tpfa_parallel-00028.pvtu
                        --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_co2_tpfa params.input -Problem.Name test_co2_tpfa_parallel")

# build target for the CO2 non-isothermal test problem
add_executable(test_co2ni_box EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_co2ni_box PUBLIC TYPETAG=HeterogeneousNIBox PUBLIC ISOTHERMAL=0 )

add_executable(test_co2ni_tpfa EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_co2ni_tpfa PUBLIC TYPETAG=HeterogeneousNICCTpfa PUBLIC ISOTHERMAL=0 )

# sequential box
dumux_add_test(NAME test_co2ni_box_sequential
               TARGET test_co2ni_box
               LABELS porousmediumflow co2 co2ni
               COMPILE_DEFINITIONS TYPETAG=HeterogeneousNIBox ISOTHERMAL=0
               CMAKE_GUARD "( dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_co2ni_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_co2ni_box-00028.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_co2ni_box paramsni.input -Problem.Name test_co2ni_box"
                        --zeroThreshold {"porosity":1})

# sequential tpfa
dumux_add_test(NAME test_co2ni_tpfa_sequential
               TARGET test_co2ni_tpfa
               LABELS porousmediumflow co2 co2ni
               CMAKE_GUARD "( dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_co2ni_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_co2ni_tpfa-00028.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_co2ni_tpfa paramsni.input -Problem.Name test_co2ni_tpfa")

# parallel box
dumux_add_test(NAME test_co2ni_box_parallel
               TARGET test_co2ni_box
               LABELS porousmediumflow co2 co2ni parallel
               TIMEOUT 1200
               CMAKE_GUARD "( dune-alugrid_FOUND AND MPI_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy --zeroThreshold {"process rank":100,"porosity":1}
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_co2ni_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_co2ni_box_parallel-00028.pvtu
                        --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_co2ni_box paramsni.input -Problem.Name test_co2ni_box_parallel")

# parallel tpfa
dumux_add_test(NAME test_co2ni_tpfa_parallel
               TARGET test_co2ni_tpfa
               LABELS porousmediumflow co2 co2ni parallel
               TIMEOUT 1200
               CMAKE_GUARD "( dune-alugrid_FOUND AND MPI_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy --zeroThreshold {"process rank":100}
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_co2ni_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_co2ni_tpfa_parallel-00028.pvtu
                        --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_co2ni_tpfa paramsni.input -Problem.Name test_co2ni_tpfa_parallel")
