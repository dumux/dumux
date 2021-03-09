add_executable(test_md_embedded_2d3d_fracture1p_tpfa EXCLUDE_FROM_ALL main.cc)

dumux_add_test(NAME test_md_embedded_2d3d_fracture1p_tpfa_linearsolver
              LABELS multidomain multidomain_embedded 1p
              TARGET test_md_embedded_2d3d_fracture1p_tpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_2d3d_fracture1p_tpfa_2d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_2d3d_fracture1p_tpfa_linearsolver_2d-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_2d3d_fracture1p_tpfa_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_2d3d_fracture1p_tpfa_linearsolver_3d-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_2d3d_fracture1p_tpfa params.input \
                                   -Vtk.OutputName test_md_embedded_2d3d_fracture1p_tpfa_linearsolver -Problem.SolverType linear")

dumux_add_test(NAME test_md_embedded_2d3d_fracture1p_tpfa_nonlinearsolver
              LABELS multidomain multidomain_embedded 1p
              TARGET test_md_embedded_2d3d_fracture1p_tpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_2d3d_fracture1p_tpfa_2d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_2d3d_fracture1p_tpfa_nonlinearsolver_2d-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_2d3d_fracture1p_tpfa_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_2d3d_fracture1p_tpfa_nonlinearsolver_3d-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_2d3d_fracture1p_tpfa params.input \
                                   -Vtk.OutputName test_md_embedded_2d3d_fracture1p_tpfa_nonlinearsolver -Problem.SolverType nonlinear")
dune_symlink_to_source_files(FILES "fracture.msh" "params.input")
