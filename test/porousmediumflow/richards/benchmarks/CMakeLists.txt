dune_symlink_to_source_files(FILES
params_evaporation_sand.input
params_evaporation_loam1.input
params_evaporation_loam2.input
params_evaporation_clay.input
params_infiltration_sand.input
params_infiltration_loam.input
params_infiltration_clay.input
run_and_plot_m21.py
run_and_plot_m22.py
)

add_executable(test_richards_benchmark_tpfa EXCLUDE_FROM_ALL main.cc)

dumux_add_test(NAME test_richards_benchmark_infiltration_tpfa
               TARGET test_richards_benchmark_tpfa
               CMAKE_GUARD HAVE_UMFPACK
               LABELS porousmediumflow richards
               TIMEOUT 3600
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzyData
                        --relative 0.05
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_richards_benchmark_tpfa_theta_deltaeta_ana_loam_0.dat
                                ${CMAKE_CURRENT_BINARY_DIR}/theta_deltaeta_ana_loam_0.dat
                                ${CMAKE_SOURCE_DIR}/test/references/test_richards_benchmark_tpfa_theta_deltaeta_num_loam_1.dat
                                ${CMAKE_CURRENT_BINARY_DIR}/theta_deltaeta_num_loam_1.dat
                                ${CMAKE_SOURCE_DIR}/test/references/test_richards_benchmark_tpfa_theta_deltaeta_num_loam_2.dat
                                ${CMAKE_CURRENT_BINARY_DIR}/theta_deltaeta_num_loam_2.dat
                                ${CMAKE_SOURCE_DIR}/test/references/test_richards_benchmark_tpfa_theta_deltaeta_num_loam_3.dat
                                ${CMAKE_CURRENT_BINARY_DIR}/theta_deltaeta_num_loam_3.dat
                                ${CMAKE_SOURCE_DIR}/test/references/test_richards_benchmark_tpfa_theta_deltaeta_ana_clay_0.dat
                                ${CMAKE_CURRENT_BINARY_DIR}/theta_deltaeta_ana_clay_0.dat
                                ${CMAKE_SOURCE_DIR}/test/references/test_richards_benchmark_tpfa_theta_deltaeta_num_clay_1.dat
                                ${CMAKE_CURRENT_BINARY_DIR}/theta_deltaeta_num_clay_1.dat
                                ${CMAKE_SOURCE_DIR}/test/references/test_richards_benchmark_tpfa_theta_deltaeta_num_clay_2.dat
                                ${CMAKE_CURRENT_BINARY_DIR}/theta_deltaeta_num_clay_2.dat
                                ${CMAKE_SOURCE_DIR}/test/references/test_richards_benchmark_tpfa_theta_deltaeta_num_clay_3.dat
                                ${CMAKE_CURRENT_BINARY_DIR}/theta_deltaeta_num_clay_3.dat
                                ${CMAKE_SOURCE_DIR}/test/references/test_richards_benchmark_tpfa_theta_deltaeta_ana_sand_0.dat
                                ${CMAKE_CURRENT_BINARY_DIR}/theta_deltaeta_ana_sand_0.dat
                                ${CMAKE_SOURCE_DIR}/test/references/test_richards_benchmark_tpfa_theta_deltaeta_num_sand_1.dat
                                ${CMAKE_CURRENT_BINARY_DIR}/theta_deltaeta_num_sand_1.dat
                                ${CMAKE_SOURCE_DIR}/test/references/test_richards_benchmark_tpfa_theta_deltaeta_num_sand_2.dat
                                ${CMAKE_CURRENT_BINARY_DIR}/theta_deltaeta_num_sand_2.dat
                                ${CMAKE_SOURCE_DIR}/test/references/test_richards_benchmark_tpfa_theta_deltaeta_num_sand_3.dat
                                ${CMAKE_CURRENT_BINARY_DIR}/theta_deltaeta_num_sand_3.dat

                        --command "${CMAKE_CURRENT_BINARY_DIR}/run_and_plot_m21.py")

dumux_add_test(NAME test_richards_benchmark_evaporation_tpfa
               TARGET test_richards_benchmark_tpfa
               CMAKE_GUARD HAVE_UMFPACK
               LABELS porousmediumflow richards
               TIMEOUT 3600
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzyData
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_richards_benchmark_tpfa_rate_analytical_clay.dat
                                ${CMAKE_CURRENT_BINARY_DIR}/rate_analytical_clay.dat
                                ${CMAKE_SOURCE_DIR}/test/references/test_richards_benchmark_tpfa_rate_analytical_loam1.dat
                                ${CMAKE_CURRENT_BINARY_DIR}/rate_analytical_loam1.dat
                                ${CMAKE_SOURCE_DIR}/test/references/test_richards_benchmark_tpfa_rate_analytical_loam2.dat
                                ${CMAKE_CURRENT_BINARY_DIR}/rate_analytical_loam2.dat
                                ${CMAKE_SOURCE_DIR}/test/references/test_richards_benchmark_tpfa_rate_analytical_sand.dat
                                ${CMAKE_CURRENT_BINARY_DIR}/rate_analytical_sand.dat
                                ${CMAKE_SOURCE_DIR}/test/references/test_richards_benchmark_tpfa_rate_actual_clay.dat
                                ${CMAKE_CURRENT_BINARY_DIR}/rate_actual_clay.dat
                                ${CMAKE_SOURCE_DIR}/test/references/test_richards_benchmark_tpfa_rate_actual_loam1.dat
                                ${CMAKE_CURRENT_BINARY_DIR}/rate_actual_loam1.dat
                                ${CMAKE_SOURCE_DIR}/test/references/test_richards_benchmark_tpfa_rate_actual_loam2.dat
                                ${CMAKE_CURRENT_BINARY_DIR}/rate_actual_loam2.dat
                                ${CMAKE_SOURCE_DIR}/test/references/test_richards_benchmark_tpfa_rate_actual_sand.dat
                                ${CMAKE_CURRENT_BINARY_DIR}/rate_actual_sand.dat
                        --command "${CMAKE_CURRENT_BINARY_DIR}/run_and_plot_m22.py")
