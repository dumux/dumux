add_input_file_links()

dune_add_test(NAME test_stokes1p2cdarcy2p2chorizontal
              SOURCES test_stokes1p2cdarcy2p2chorizontal.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_stokes1p2cdarcy2p2chorizontal_stokes-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_stokes1p2cdarcy2p2chorizontal_stokes-00040.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_stokes1p2cdarcy2p2chorizontal_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_stokes1p2cdarcy2p2chorizontal_darcy-00040.vtu

                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes1p2cdarcy2p2chorizontal test_stokes1p2cdarcy2p2chorizontal.input")

dune_add_test(NAME test_stokes1p2cnidarcy2p2cnihorizontal
              SOURCES test_stokes1p2cdarcy2p2chorizontal.cc
              COMPILE_DEFINITIONS NONISOTHERMAL=1
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_stokes1p2cnidarcy2p2cnihorizontal_stokes-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_stokes1p2cnidarcy2p2cnihorizontal_stokes-00039.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_stokes1p2cnidarcy2p2cnihorizontal_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_stokes1p2cnidarcy2p2cnihorizontal_darcy-00039.vtu

                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes1p2cnidarcy2p2cnihorizontal test_stokes1p2cdarcy2p2chorizontal.input -Problem.Name test_stokes1p2cnidarcy2p2cnihorizontal")
