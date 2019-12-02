add_input_file_links()

dumux_add_test(NAME test_md_boundary_darcy1p2c_stokes1p2c_maxwellstefan
              LABELS multidomain multidomain_boundary freeflow 1pnc
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS DIFFUSIONTYPE=MaxwellStefansLaw<TypeTag>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p2c_stokes1p2c_comparison_stokes-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2c_stokes1p2c_comparison_ms_stokes-00020.vtu
                                    ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p2c_stokes1p2c_comparison_darcy-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2c_stokes1p2c_comparison_ms_darcy-00020.vtu

                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2c_stokes1p2c_maxwellstefan params.input
                                        -Problem.Name test_md_boundary_darcy1p2c_stokes1p2c_comparison_ms")

dumux_add_test(NAME test_md_boundary_darcy1p2c_stokes1p2c_fickslaw
              LABELS multidomain multidomain_boundary freeflow 1pnc
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS DIFFUSIONTYPE=FicksLaw<TypeTag>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p2c_stokes1p2c_comparison_stokes-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2c_stokes1p2c_comparison_f_stokes-00020.vtu
                                    ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p2c_stokes1p2c_comparison_darcy-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2c_stokes1p2c_comparison_f_darcy-00020.vtu

                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p2c_stokes1p2c_fickslaw params.input
                                        -Problem.Name test_md_boundary_darcy1p2c_stokes1p2c_comparison_f")
