add_input_file_links()

dumux_add_test(NAME test_md_boundary_darcy2p2c_stokes1p2c_horizontal
              LABELS multidomain multidomain_boundary stokesdarcy 2p2c navierstokes
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2c_stokes1p2c_horizontal_stokes-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2c_stokes1p2c_horizontal_stokes-00039.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2c_stokes1p2c_horizontal_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2c_stokes1p2c_horizontal_darcy-00039.vtu

                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2c_stokes1p2c_horizontal params.input
                             -Vtk.OutputName test_md_boundary_darcy2p2c_stokes1p2c_horizontal")

dumux_add_test(NAME test_md_boundary_darcy2p2cni_stokes1p2cni_horizontal
              LABELS multidomain multidomain_boundary stokesdarcy 2p2c navierstokes
              SOURCES main.cc
              COMPILE_DEFINITIONS NONISOTHERMAL=1
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2cni_stokes1p2cni_horizontal_stokes-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_stokes1p2cni_horizontal_stokes-00039.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2cni_stokes1p2cni_horizontal_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_stokes1p2cni_horizontal_darcy-00039.vtu

                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_stokes1p2cni_horizontal params_nonisothermal.input
                             -Vtk.OutputName test_md_boundary_darcy2p2cni_stokes1p2cni_horizontal")
