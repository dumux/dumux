dumux_add_test(SOURCES test_thermalconductivity.cc
              LABELS unit material
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzyData --delimiter " "
                      --files ${CMAKE_SOURCE_DIR}/test/references/thermalconductivitysomerton-reference.dat
                              ${CMAKE_CURRENT_BINARY_DIR}/somerton_lambda_eff.dat
                              ${CMAKE_SOURCE_DIR}/test/references/thermalconductivityjohansen-reference.dat
                              ${CMAKE_CURRENT_BINARY_DIR}/johansen_lambda_eff.dat
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_thermalconductivity")

dumux_add_test(SOURCES test_material_2p_vangenuchten.cc
               LABELS unit material
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzyData --delimiter " "
                        --files  ${CMAKE_SOURCE_DIR}/test/references/test_pcsw_vangenuchten.dat
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_pcsw_vangenuchten.dat
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_material_2p_vangenuchten")

dumux_add_test(SOURCES test_material_2p_brookscorey.cc
               LABELS unit material
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzyData --delimiter " "
                        --files  ${CMAKE_SOURCE_DIR}/test/references/test_pcsw_brookscorey.dat
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_pcsw_brookscorey.dat
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_material_2p_brookscorey")


dune_add_test(SOURCES test_material_2p_spline.cc
              LABELS unit material)

dune_add_test(SOURCES test_material_2p_dataspline.cc
              LABELS unit material)

dune_symlink_to_source_files(FILES test_material_2p_spline.input test_material_2p_dataspline.input)
