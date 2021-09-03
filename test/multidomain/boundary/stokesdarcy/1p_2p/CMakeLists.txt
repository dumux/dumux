add_input_file_links()

add_executable(test_md_boundary_pm2ptpfa_ff1p_vertical EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_md_boundary_pm2ptpfa_ff1p_vertical PUBLIC DARCYTYPETAG=DarcyTwoPCC)

add_executable(test_md_boundary_pm2pbox_ff1p_vertical EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_md_boundary_pm2pbox_ff1p_vertical PUBLIC DARCYTYPETAG=DarcyTwoPBox)

dumux_add_test(NAME test_md_boundary_darcy2ptpfa_stokes1p_vertical
              TARGET test_md_boundary_pm2ptpfa_ff1p_vertical
              LABELS multidomain multidomain_boundary stokesdarcy 1p navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2ptpfa_stokes1p_vertical_stokes-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2ptpfa_stokes1p_vertical_stokes-00024.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2ptpfa_stokes1p_vertical_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2ptpfa_stokes1p_vertical_darcy-00024.vtu

                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_pm2ptpfa_ff1p_vertical params.input
                             -Vtk.OutputName test_md_boundary_darcy2ptpfa_stokes1p_vertical"
                             --zeroThreshold {"velocity_liq \(m/s\)_0":1e-14,"velocity_gas \(m/s\)":1e-15})

dumux_add_test(NAME test_md_boundary_darcy2pbox_stokes1p_vertical
              TARGET test_md_boundary_pm2pbox_ff1p_vertical
              LABELS multidomain multidomain_boundary stokesdarcy 1p navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                              --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2pbox_stokes1p_vertical_stokes-reference.vtu
                                      ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2pbox_stokes1p_vertical_stokes-00024.vtu
                                      ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2pbox_stokes1p_vertical_darcy-reference.vtu
                                      ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2pbox_stokes1p_vertical_darcy-00024.vtu
                              --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_pm2pbox_ff1p_vertical params.input
                              -Vtk.OutputName test_md_boundary_darcy2pbox_stokes1p_vertical"
                              --zeroThreshold {"velocity_liq \(m/s\)_0":2e-14,"velocity_gas \(m/s\)":5e-13})
