add_input_file_links()

dune_add_test(NAME test_stokes1p2cdarcy1p2chorizontal
             SOURCES test_stokes1p2cdarcy1p2chorizontal.cc
             CMAKE_GUARD HAVE_UMFPACK
             COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
             CMD_ARGS       --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_stokes1p2cdarcy1p2chorizontal_stokes-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/test_stokes1p2cdarcy1p2chorizontal_stokes-00020.vtu
                                    ${CMAKE_SOURCE_DIR}/test/references/test_stokes1p2cdarcy1p2chorizontal_darcy-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/test_stokes1p2cdarcy1p2chorizontal_darcy-00020.vtu

                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes1p2cdarcy1p2chorizontal test_stokes1p2cdarcy1p2chorizontal.input")
