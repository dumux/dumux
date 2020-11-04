add_subdirectory(convergencetest)

add_input_file_links()

add_executable(test_md_boundary_darcy1p_stokes1p EXCLUDE_FROM_ALL main.cc)

dumux_add_test(NAME test_md_boundary_darcy1p_stokes1p_horizontal
              LABELS multidomain multidomain_boundary stokesdarcy 1p navierstokes
              TARGET test_md_boundary_darcy1p_stokes1p
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_stokes1p_horizontal_stokes-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_stokes1p_horizontal_stokes-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_stokes1p_horizontal_darcy-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_stokes1p_horizontal_darcy-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_stokes1p params_horizontalflow.input
                                   -Vtk.OutputName test_md_boundary_darcy1p_stokes1p_horizontal")

dumux_add_test(NAME test_md_boundary_darcy1p_stokes1p_vertical
              LABELS multidomain multidomain_boundary stokesdarcy 1p navierstokes
              TARGET test_md_boundary_darcy1p_stokes1p
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_stokes1p_vertical_stokes-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_stokes1p_vertical_stokes-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_stokes1p_vertical_darcy-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_stokes1p_vertical_darcy-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_stokes1p params_verticalflow.input
                                   -Vtk.OutputName test_md_boundary_darcy1p_stokes1p_vertical"
                        --zeroThreshold {"velocity_liq \(m/s\)_0":6e-17})
