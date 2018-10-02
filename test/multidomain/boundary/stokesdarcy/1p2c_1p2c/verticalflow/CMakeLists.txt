add_input_file_links()

add_executable(test_stokes1p2cdarcy1p2cvertical EXCLUDE_FROM_ALL test_stokes1p2cdarcy1p2cvertical.cc)

dune_add_test(NAME test_stokes1p2cdarcy1p2cvertical_diffusion
              TARGET test_stokes1p2cdarcy1p2cvertical
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --zeroThreshold {"velocity_liq \(m/s\)":1e-20}
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_stokes1p2cdarcy1p2cvertical_diffusion_stokes-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_stokes1p2cdarcy1p2cvertical_diffusion_stokes-00003.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_stokes1p2cdarcy1p2cvertical_diffusion_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_stokes1p2cdarcy1p2cvertical_diffusion_darcy-00003.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes1p2cdarcy1p2cvertical test_stokes1p2cdarcy1p2cvertical_diffusion.input")

dune_add_test(NAME test_stokes1p2cdarcy1p2cvertical_advection
              TARGET test_stokes1p2cdarcy1p2cvertical
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --zeroThreshold {"velocity_liq \(m/s\)":1e-15}
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_stokes1p2cdarcy1p2cvertical_advection_stokes-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_stokes1p2cdarcy1p2cvertical_stokes-00030.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_stokes1p2cdarcy1p2cvertical_advection_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_stokes1p2cdarcy1p2cvertical_darcy-00030.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes1p2cdarcy1p2cvertical test_stokes1p2cdarcy1p2cvertical.input")
