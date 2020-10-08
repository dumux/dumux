add_input_file_links()

execute_process(COMMAND ${CMAKE_COMMAND} "-E" "create_symlink" "./../convergencetest/convergencetest.py" "${CMAKE_CURRENT_BINARY_DIR}/convergencetest.py")

dune_add_test(NAME test_md_boundary_darcy1p_stokes1p_bj
              SOURCES main.cc
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND convergencetest.py
              CMD_ARGS test_md_boundary_darcy1p_stokes1p_bj params_bj.input)

dune_add_test(NAME test_md_boundary_darcy1p_stokes1p_newic
              SOURCES main.cc
              COMPILE_DEFINITIONS NEWIC=1
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND convergencetest.py
              CMD_ARGS test_md_boundary_darcy1p_stokes1p_newic params_newic.input)