# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()
dune_symlink_to_source_files(FILES grids)

# build target for the CO2 test problem
# Ignore the porosity for all box models since it is defined element-wise in these tests
# but the default 2p2c implementation outputs porosity per vertex.
# Depending on the order of the elements, the porosity would differ in these cases.

dumux_add_test(NAME test_co2_box
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=HeterogeneousBox
               LABELS porousmediumflow co2
               TIMEOUT 3600
               CMAKE_GUARD "( dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_co2_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_co2_box-00028.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_co2_box params.input -Problem.Name test_co2_box"
                        --ignore "porosity")

dumux_add_test(NAME test_co2_tpfa
               LABELS porousmediumflow co2
               TIMEOUT 3600
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=HeterogeneousCCTpfa
               CMAKE_GUARD "( dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_co2_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_co2_tpfa-00028.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_co2_tpfa params.input -Problem.Name test_co2_tpfa")

dumux_add_test(NAME test_co2_mpfa
               LABELS porousmediumflow co2
               TIMEOUT 3600
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=HeterogeneousCCMpfa
               CMAKE_GUARD "( dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_co2_mpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_co2_mpfa-00028.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_co2_mpfa params.input -Problem.Name test_co2_mpfa")

dumux_add_test(NAME test_co2_mpfa_fluxvarscache
               LABELS porousmediumflow co2
               TIMEOUT 3600
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=HeterogeneousCCMpfa
               COMPILE_DEFINITIONS FLUXVARSCACHE=1
               CMAKE_GUARD "( dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_co2_mpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_co2_mpfa_fluxvarscache-00028.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_co2_mpfa_fluxvarscache params.input -Problem.Name test_co2_mpfa_fluxvarscache")


dumux_add_test(NAME test_co2_box_parallel
               TARGET test_co2_box
               LABELS porousmediumflow co2 parallel
               TIMEOUT 3600
               CMAKE_GUARD "( dune-alugrid_FOUND AND MPI_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy --ignore "process rank" "porosity"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_co2_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_co2_box_parallel-00028.pvtu
                        --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_co2_box params.input -Problem.Name test_co2_box_parallel")

dumux_add_test(NAME test_co2_tpfa_parallel
               TARGET test_co2_tpfa
               LABELS porousmediumflow co2 parallel
               TIMEOUT 3600
               CMAKE_GUARD "( dune-alugrid_FOUND AND MPI_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy --ignore "process rank"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_co2_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_co2_tpfa_parallel-00028.pvtu
                        --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_co2_tpfa params.input -Problem.Name test_co2_tpfa_parallel")


dumux_add_test(NAME test_co2ni_box
               SOURCES main.cc
               LABELS porousmediumflow co2 co2ni
               TIMEOUT 3600
               COMPILE_DEFINITIONS TYPETAG=HeterogeneousNIBox ISOTHERMAL=0
               CMAKE_GUARD "( dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_co2ni_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_co2ni_box-00028.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_co2ni_box paramsni.input -Problem.Name test_co2ni_box"
                        --ignore "porosity")

dumux_add_test(NAME test_co2ni_tpfa
               SOURCES main.cc
               LABELS porousmediumflow co2 co2ni
               TIMEOUT 3600
               COMPILE_DEFINITIONS TYPETAG=HeterogeneousNICCTpfa PUBLIC ISOTHERMAL=0
               CMAKE_GUARD "( dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_co2ni_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_co2ni_tpfa-00028.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_co2ni_tpfa paramsni.input -Problem.Name test_co2ni_tpfa")

dumux_add_test(NAME test_co2ni_box_parallel
               TARGET test_co2ni_box
               LABELS porousmediumflow co2 co2ni parallel
               TIMEOUT 3600
               CMAKE_GUARD "( dune-alugrid_FOUND AND MPI_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy --ignore "process rank" "porosity"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_co2ni_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_co2ni_box_parallel-00028.pvtu
                        --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_co2ni_box paramsni.input -Problem.Name test_co2ni_box_parallel")

dumux_add_test(NAME test_co2ni_tpfa_parallel
               TARGET test_co2ni_tpfa
               LABELS porousmediumflow co2 co2ni parallel
               TIMEOUT 3600
               CMAKE_GUARD "( dune-alugrid_FOUND AND MPI_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy --ignore "process rank"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_co2ni_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_co2ni_tpfa_parallel-00028.pvtu
                        --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_co2ni_tpfa paramsni.input -Problem.Name test_co2ni_tpfa_parallel")
