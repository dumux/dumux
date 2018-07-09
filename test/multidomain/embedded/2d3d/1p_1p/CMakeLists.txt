add_executable(test_embedded_fracture1p_tpfa EXCLUDE_FROM_ALL test_1p_1p.cc)

dune_add_test(NAME test_embedded_fracture1p_tpfa_linearsolver
              TARGET test_embedded_fracture1p_tpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_embedded_fracture_tpfa_2d.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_embedded_fracture_tpfa_linearsolver_2d-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_embedded_fracture_tpfa_3d.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_embedded_fracture_tpfa_linearsolver_3d-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_embedded_fracture1p_tpfa test_embedded_fracture1p.input \
                                   -Problem.Name test_embedded_fracture_tpfa_linearsolver -Problem.SolverType linear")

dune_add_test(NAME test_embedded_fracture1p_tpfa_nonlinearsolver
              TARGET test_embedded_fracture1p_tpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_embedded_fracture_tpfa_2d.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_embedded_fracture_tpfa_nonlinearsolver_2d-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_embedded_fracture_tpfa_3d.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_embedded_fracture_tpfa_nonlinearsolver_3d-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_embedded_fracture1p_tpfa test_embedded_fracture1p.input \
                                   -Problem.Name test_embedded_fracture_tpfa_nonlinearsolver -Problem.SolverType nonlinear")
dune_symlink_to_source_files(FILES "fracture.msh" "test_embedded_fracture1p.input")
set(CMAKE_BUILD_TYPE Release)
