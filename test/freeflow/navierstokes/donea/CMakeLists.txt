add_input_file_links()
dune_symlink_to_source_files(FILES grids)

dumux_add_test(NAME test_ff_stokes_donea_momentum
               LABELS freeflow navierstokes
               SOURCES main_momentum.cc
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/donea_momentum_1.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum")

dumux_add_test(NAME test_ff_stokes_donea_momentum_parallel
               LABELS freeflow navierstokes
               TARGET test_ff_stokes_donea_momentum
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI AND ( DUNE_ISTL_VERSION VERSION_GREATER_EQUAL 2.8 ) )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy --zeroThreshold {"rank":100}
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0004-donea_momentum_1.pvtu
                             --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum")

dumux_add_test(NAME test_donea_momentum_diamond_rt
               SOURCES main_momentum_diamond.cc)

dumux_add_test(NAME test_donea_momentum_diamond_cr
               SOURCES main_momentum_diamond.cc
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>)

dumux_add_test(NAME test_ff_stokes_donea_nocaching
              LABELS freeflow navierstokes
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS ENABLECACHING=0
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_nocaching-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_nocaching params.input
                             -Problem.Name test_ff_stokes_donea_nocaching -Problem.UseNeumann false")

dumux_add_test(NAME test_ff_stokes_donea_donut_nocaching
              LABELS freeflow navierstokes
              TARGET test_ff_stokes_donea_nocaching
              CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_donut-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_donut_nocaching-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_nocaching params.input
                             -Grid.File grids/donut.dgf
                             -Problem.Name test_ff_stokes_donea_donut_nocaching -Problem.UseNeumann false")

dumux_add_test(NAME test_ff_stokes_donea_donut_twisted_nocaching
              LABELS freeflow navierstokes
              TARGET test_ff_stokes_donea_nocaching
              CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_donut-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_donut_twisted_nocaching-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_nocaching params.input
                             -Grid.File grids/donut_twisted.dgf
                             -Problem.Name test_ff_stokes_donea_donut_twisted_nocaching -Problem.UseNeumann false"
                             --zeroThreshold {"connectivity":1e6})

dumux_add_test(NAME test_ff_stokes_donea
              LABELS freeflow navierstokes
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS ENABLECACHING=1
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea params.input
                             -Problem.Name test_ff_stokes_donea -Problem.UseNeumann false")

dumux_add_test(NAME test_ff_stokes_donea_donut
               LABELS freeflow navierstokes
               TARGET test_ff_stokes_donea
               CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS       --script fuzzy
                               --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_donut-reference.vtu
                                       ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_donut-00001.vtu
                               --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea params.input
                               -Grid.File grids/donut.dgf
                               -Problem.Name test_ff_stokes_donea_donut -Problem.UseNeumann false")

dumux_add_test(NAME test_ff_stokes_donea_donut_twisted
               LABELS freeflow navierstokes
               TARGET test_ff_stokes_donea
               CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS       --script fuzzy
                               --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_donut-reference.vtu
                                       ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_donut_twisted-00001.vtu
                               --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea params.input
                               -Grid.File grids/donut_twisted.dgf
                               -Problem.Name test_ff_stokes_donea_donut_twisted -Problem.UseNeumann false"
                               --zeroThreshold {"connectivity":1e6})
