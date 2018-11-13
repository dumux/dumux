add_input_file_links()

dune_add_test(NAME test_md_boundary_darcy1p_stokes1p_horizontal
              LABELS multidomain
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_stokes1p_horizontal_stokes-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_stokes1p_horizontal_stokes-00001.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_stokes1p_horizontal_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_stokes1p_horizontal_darcy-00001.vtu

                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_stokes1p_horizontal params.input
                             -Vtk.OutputName test_md_boundary_darcy1p_stokes1p_horizontal")
