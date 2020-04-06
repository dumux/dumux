dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_md_poromechanics_el1p
               LABELS multidomain poromechanics 1p poroelastic
               SOURCES main.cc
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_poromechanics_el1p_1p-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_poromechanics_el1p_onep-00001.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_poromechanics_el1p_poroelastic-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_md_poromechanics_el1p_poroelastic-00001.vtu
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_poromechanics_el1p params.input
                                                               -Vtk.OutputName test_md_poromechanics_el1p"
                         --zeroThreshold {"u":1e-14})

dumux_add_test(NAME test_md_poromechanics_el1p_parallel
               TARGET test_md_poromechanics_el1p
               LABELS multidomain poromechanics 1p poroelastic
               CMAKE_GUARD MPI_FOUND
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_poromechanics_el1p_1p-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_md_poromechanics_el1p_parallel_onep-00001.pvtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_poromechanics_el1p_poroelastic-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_md_poromechanics_el1p_parallel_poroelastic-00001.pvtu
                         --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_md_poromechanics_el1p params.input
                                    -Grid.Overlap 1 -Vtk.OutputName test_md_poromechanics_el1p_parallel"
                         --zeroThreshold {"u":1e-14,"process rank":100})
