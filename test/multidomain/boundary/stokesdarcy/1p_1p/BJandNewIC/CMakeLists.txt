add_input_file_links()

add_executable(test_md_boundary_darcy1p_stokes1p_New EXCLUDE_FROM_ALL main.cc)

dumux_add_test(NAME test_md_boundary_darcy1p_stokes1p_BJ
              LABELS multidomain multidomain_boundary stokesdarcy 1p navierstokes
              TARGET test_md_boundary_darcy1p_stokes1p_New
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_stokes1p_BJ_stokes-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_stokes1p_BJ_stokes-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_stokes1p_BJ_darcy-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_stokes1p_BJ_darcy-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_stokes1p_New params_BJ.input
                                   -Vtk.OutputName test_md_boundary_darcy1p_stokes1p_BJ")

dumux_add_test(NAME test_md_boundary_darcy1p_stokes1p_NewIC
              LABELS multidomain multidomain_boundary stokesdarcy 1p navierstokes
              TARGET test_md_boundary_darcy1p_stokes1p_New
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_stokes1p_NewIC_stokes-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_stokes1p_NewIC_stokes-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_stokes1p_NewIC_darcy-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_stokes1p_NewIC_darcy-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_stokes1p_New params_newIC.input
                                   -Vtk.OutputName test_md_boundary_darcy1p_stokes1p_NewIC"
                        --zeroThreshold {"velocity_liq \(m/s\)_0":6e-17})
