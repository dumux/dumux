add_input_file_links()

dune_add_test(NAME test_stokes1pdarcy2pvertical
              SOURCES test_stokes1pdarcy2pvertical.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_stokes1pdarcy2pvertical_stokes-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_stokes1pdarcy2pvertical_stokes-00024.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_stokes1pdarcy2pvertical_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_stokes1pdarcy2pvertical_darcy-00024.vtu

                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes1pdarcy2pvertical test_stokes1pdarcy2pvertical.input"
                             --zeroThreshold {"velocity_liq \(m/s\)_0":1e-14,"velocity_gas \(m/s\)":1e-15})
