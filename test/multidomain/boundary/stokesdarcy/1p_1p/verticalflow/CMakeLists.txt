add_input_file_links()

dune_add_test(NAME test_stokes1pdarcy1pvertical
              SOURCES test_stokes1pdarcy1pvertical.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_stokes1pdarcy1pvertical_stokes-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_stokes1pdarcy1pvertical_stokes-00001.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_stokes1pdarcy1pvertical_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_stokes1pdarcy1pvertical_darcy-00001.vtu

                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes1pdarcy1pvertical test_stokes1pdarcy1pvertical.input"
                             --zeroThreshold {"velocity_liq \(m/s\)_0":6e-17})
