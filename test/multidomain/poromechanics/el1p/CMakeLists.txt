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

dumux_add_test(NAME test_md_poromechanics_el1p_parallel_yasp
               TARGET test_md_poromechanics_el1p
               LABELS multidomain poromechanics 1p poroelastic
               CMAKE_GUARD MPI_FOUND
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_poromechanics_el1p_1p-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_md_poromechanics_el1p_parallel_yasp_onep-00001.pvtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_poromechanics_el1p_poroelastic-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_md_poromechanics_el1p_parallel_yasp_poroelastic-00001.pvtu
                         --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_md_poromechanics_el1p params.input
                                    -Grid.Overlap 1 -Vtk.OutputName test_md_poromechanics_el1p_parallel_yasp"
                         --zeroThreshold {"u":1e-14,"process rank":100})

dumux_add_test(NAME test_md_poromechanics_el1p_parallel_ug
               SOURCES main.cc
               LABELS multidomain poromechanics 1p poroelastic
               CMAKE_GUARD "( MPI_FOUND AND HAVE_UG )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::UGGrid<2>
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_poromechanics_el1p_1p-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_md_poromechanics_el1p_parallel_ug_onep-00001.pvtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_poromechanics_el1p_poroelastic-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_md_poromechanics_el1p_parallel_ug_poroelastic-00001.pvtu
                         --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_md_poromechanics_el1p_parallel_ug
                                    params.input -Vtk.OutputName test_md_poromechanics_el1p_parallel_ug"
                         --zeroThreshold {"u":1e-14,"process rank":100})

dumux_add_test(NAME test_md_poromechanics_el1p_parallel_alu
               SOURCES main.cc
               LABELS multidomain poromechanics 1p poroelastic
               CMAKE_GUARD "( MPI_FOUND AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_md_poromechanics_el1p_1p-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_md_poromechanics_el1p_parallel_alu_onep-00001.pvtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_md_poromechanics_el1p_poroelastic-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_md_poromechanics_el1p_parallel_alu_poroelastic-00001.pvtu
                         --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_md_poromechanics_el1p_parallel_alu
                                    params.input -Vtk.OutputName test_md_poromechanics_el1p_parallel_alu"
                         --zeroThreshold {"u":1e-14,"process rank":100})
