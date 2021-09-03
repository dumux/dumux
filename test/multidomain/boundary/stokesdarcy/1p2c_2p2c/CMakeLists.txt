add_input_file_links()

add_executable(test_md_boundary_darcy2p2ctpfa_stokes1p2c_horizontal EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_md_boundary_darcy2p2ctpfa_stokes1p2c_horizontal PUBLIC DARCYTYPETAG=DarcyTwoPTwoCTpfa)

add_executable(test_md_boundary_darcy2p2cnitpfa_stokes1p2cni_horizontal EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_md_boundary_darcy2p2cnitpfa_stokes1p2cni_horizontal PUBLIC "DARCYTYPETAG=DarcyTwoPTwoCTpfa" "NONISOTHERMAL=1")

add_executable(test_md_boundary_darcy2p2cbox_stokes1p2c_horizontal EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_md_boundary_darcy2p2cbox_stokes1p2c_horizontal PUBLIC DARCYTYPETAG=DarcyTwoPTwoCBox)

add_executable(test_md_boundary_darcy2p2cnibox_stokes1p2cni_horizontal EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_md_boundary_darcy2p2cnibox_stokes1p2cni_horizontal PUBLIC "DARCYTYPETAG=DarcyTwoPTwoCBox" "NONISOTHERMAL=1")

dumux_add_test(NAME test_md_boundary_darcy2p2ctpfa_stokes1p2c_horizontal
              TARGET test_md_boundary_darcy2p2ctpfa_stokes1p2c_horizontal
              LABELS multidomain multidomain_boundary stokesdarcy 2p2c navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2ctpfa_stokes1p2c_horizontal_stokes-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2ctpfa_stokes1p2c_horizontal_stokes-00039.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2ctpfa_stokes1p2c_horizontal_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2ctpfa_stokes1p2c_horizontal_darcy-00039.vtu

                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2ctpfa_stokes1p2c_horizontal params.input
                             -Vtk.OutputName test_md_boundary_darcy2p2ctpfa_stokes1p2c_horizontal")

dumux_add_test(NAME test_md_boundary_darcy2p2cnitpfa_stokes1p2cni_horizontal
              TARGET test_md_boundary_darcy2p2cnitpfa_stokes1p2cni_horizontal
              LABELS multidomain multidomain_boundary stokesdarcy 2p2c navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2cnitpfa_stokes1p2cni_horizontal_stokes-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cnitpfa_stokes1p2cni_horizontal_stokes-00039.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2cnitpfa_stokes1p2cni_horizontal_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cnitpfa_stokes1p2cni_horizontal_darcy-00039.vtu

                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cnitpfa_stokes1p2cni_horizontal params_nonisothermal.input
                             -Vtk.OutputName test_md_boundary_darcy2p2cnitpfa_stokes1p2cni_horizontal")

dumux_add_test(NAME test_md_boundary_darcy2p2cbox_stokes1p2c_horizontal
              TARGET test_md_boundary_darcy2p2cbox_stokes1p2c_horizontal
              LABELS multidomain multidomain_boundary stokesdarcy 2p2c navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                              --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2cbox_stokes1p2c_horizontal_stokes-reference.vtu
                                      ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cbox_stokes1p2c_horizontal_stokes-00040.vtu
                                      ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2cbox_stokes1p2c_horizontal_darcy-reference.vtu
                                      ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cbox_stokes1p2c_horizontal_darcy-00040.vtu
                              --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cbox_stokes1p2c_horizontal params.input
                              -Vtk.OutputName test_md_boundary_darcy2p2cbox_stokes1p2c_horizontal")

dumux_add_test(NAME test_md_boundary_darcy2p2cnibox_stokes1p2cni_horizontal
              TARGET test_md_boundary_darcy2p2cnibox_stokes1p2cni_horizontal
              LABELS multidomain multidomain_boundary stokesdarcy 2p2c navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                              --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2cnibox_stokes1p2cni_horizontal_stokes-reference.vtu
                                      ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cnibox_stokes1p2cni_horizontal_stokes-00040.vtu
                                      ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2cnibox_stokes1p2cni_horizontal_darcy-reference.vtu
                                      ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cnibox_stokes1p2cni_horizontal_darcy-00040.vtu
                              --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cnibox_stokes1p2cni_horizontal params_nonisothermal.input
                              -Vtk.OutputName test_md_boundary_darcy2p2cnibox_stokes1p2cni_horizontal
                              -Newton.MaxRelativeShift  2e-5")
