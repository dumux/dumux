add_subdirectory(diffusionlawcomparison)

add_input_file_links()

add_executable(test_md_boundary_darcy1p2ctpfa_stokes1p2c EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_md_boundary_darcy1p2ctpfa_stokes1p2c PUBLIC DARCYTYPETAG=DarcyOnePTwoCTpfa)

add_executable(test_md_boundary_darcy1p2cbox_stokes1p2c EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_md_boundary_darcy1p2cbox_stokes1p2c PUBLIC DARCYTYPETAG=DarcyOnePTwoCBox)

dumux_add_test(NAME test_md_boundary_darcy1p2ctpfa_stokes1p2c_horizontal
              LABELS multidomain multidomain_boundary stokesdarcy 1pnc navierstokes
              TARGET test_md_boundary_darcy1p2ctpfa_stokes1p2c
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p2ctpfa_stokes1p2c_horizontal_stokes-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2ctpfa_stokes1p2c_horizontal_stokes-00020.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p2ctpfa_stokes1p2c_horizontal_darcy-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2ctpfa_stokes1p2c_horizontal_darcy-00020.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2ctpfa_stokes1p2c params_horizontalflow.input
                                  -Vtk.OutputName test_md_boundary_darcy1p2ctpfa_stokes1p2c_horizontal")

dumux_add_test(NAME test_md_boundary_darcy1p2ctpfa_stokes1p2c_vertical_diffusion
              LABELS multidomain multidomain_boundary stokesdarcy 1pnc navierstokes
              TARGET test_md_boundary_darcy1p2ctpfa_stokes1p2c
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --zeroThreshold {"velocity_liq \(m/s\)":1e-20}
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p2ctpfa_stokes1p2c_vertical_diffusion_stokes-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2ctpfa_stokes1p2c_vertical_diffusion_stokes-00003.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p2ctpfa_stokes1p2c_vertical_diffusion_darcy-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2ctpfa_stokes1p2c_vertical_diffusion_darcy-00003.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2ctpfa_stokes1p2c params_verticalflow_diffusion.input
                                  -Vtk.OutputName test_md_boundary_darcy1p2ctpfa_stokes1p2c_vertical_diffusion")

dumux_add_test(NAME test_md_boundary_darcy1p2ctpfa_stokes1p2c_vertical_advection
              LABELS multidomain multidomain_boundary stokesdarcy 1pnc navierstokes
              TARGET test_md_boundary_darcy1p2ctpfa_stokes1p2c
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --zeroThreshold {"velocity_liq \(m/s\)":1e-15}
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p2ctpfa_stokes1p2c_vertical_advection_stokes-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2ctpfa_stokes1p2c_vertical_advection_stokes-00030.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p2ctpfa_stokes1p2c_vertical_advection_darcy-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2ctpfa_stokes1p2c_vertical_advection_darcy-00030.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2ctpfa_stokes1p2c params_verticalflow.input
                                  -Vtk.OutputName test_md_boundary_darcy1p2ctpfa_stokes1p2c_vertical_advection")

dumux_add_test(NAME test_md_boundary_darcy1p2cbox_stokes1p2c_horizontal
              LABELS multidomain multidomain_boundary stokesdarcy 1pnc navierstokes
              TARGET test_md_boundary_darcy1p2cbox_stokes1p2c
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p2cbox_stokes1p2c_horizontal_stokes-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2cbox_stokes1p2c_horizontal_stokes-00020.vtu
                              ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p2cbox_stokes1p2c_horizontal_darcy-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2cbox_stokes1p2c_horizontal_darcy-00020.vtu
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2cbox_stokes1p2c params_horizontalflow.input
                              -Vtk.OutputName test_md_boundary_darcy1p2cbox_stokes1p2c_horizontal")

dumux_add_test(NAME test_md_boundary_darcy1p2cbox_stokes1p2c_vertical_diffusion
              LABELS multidomain multidomain_boundary stokesdarcy 1pnc navierstokes
              TARGET test_md_boundary_darcy1p2cbox_stokes1p2c
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                      --zeroThreshold {"velocity_liq \(m/s\)":1e-20}
                      --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p2cbox_stokes1p2c_vertical_diffusion_stokes-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2cbox_stokes1p2c_vertical_diffusion_stokes-00003.vtu
                              ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p2cbox_stokes1p2c_vertical_diffusion_darcy-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2cbox_stokes1p2c_vertical_diffusion_darcy-00003.vtu
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2cbox_stokes1p2c params_verticalflow_diffusion.input
                              -Vtk.OutputName test_md_boundary_darcy1p2cbox_stokes1p2c_vertical_diffusion")

dumux_add_test(NAME test_md_boundary_darcy1p2cbox_stokes1p2c_vertical_advection
              LABELS multidomain multidomain_boundary stokesdarcy 1pnc navierstokes
              TARGET test_md_boundary_darcy1p2cbox_stokes1p2c
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                      --zeroThreshold {"velocity_liq \(m/s\)":1e-15}
                      --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p2cbox_stokes1p2c_vertical_advection_stokes-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2cbox_stokes1p2c_vertical_advection_stokes-00030.vtu
                              ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p2cbox_stokes1p2c_vertical_advection_darcy-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2cbox_stokes1p2c_vertical_advection_darcy-00030.vtu
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2cbox_stokes1p2c params_verticalflow.input
                              -Vtk.OutputName test_md_boundary_darcy1p2cbox_stokes1p2c_vertical_advection")