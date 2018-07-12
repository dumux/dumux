add_input_file_links()

dune_add_test(NAME test_stokes1pdarcy1phorizontal
              SOURCES test_stokes1pdarcy1phorizontal.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_stokes1pdarcy1phorizontal_stokes-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_stokes1pdarcy1phorizontal_stokes-00001.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_stokes1pdarcy1phorizontal_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_stokes1pdarcy1phorizontal_darcy-00001.vtu

                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes1pdarcy1phorizontal test_stokes1pdarcy1phorizontal.input")
