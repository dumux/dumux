add_input_file_links()

add_executable(test_stokes1pdarcy2pvertical EXCLUDE_FROM_ALL test_stokes1pdarcy2pvertical.cc)

dune_add_test(NAME test_stokes1pdarcy2pvertical
              TARGET test_stokes1pdarcy2pvertical
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_stokes1pdarcy2pvertical_stokes-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_stokes1pdarcy2pvertical_stokes-00024.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_stokes1pdarcy2pvertical_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_stokes1pdarcy2pvertical_darcy-00024.vtu

                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes1pdarcy2pvertical test_stokes1pdarcy2pvertical.input"
                             --zeroThreshold {"velocity_SimpleH2O \(m/s\)_0":1e-14,"velocity_Air \(m/s\)":1e-15})
