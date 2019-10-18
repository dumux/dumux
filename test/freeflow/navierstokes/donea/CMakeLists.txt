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
                             -Problem.Name test_ff_stokes_donea_nocaching")

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
                             -Problem.Name test_ff_stokes_donea")

dune_symlink_to_source_files(FILES "params.input")
