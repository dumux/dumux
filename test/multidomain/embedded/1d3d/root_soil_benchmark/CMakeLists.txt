# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# link the mesh and the parameter file
dune_symlink_to_source_files(FILES params.input dap8_rootsys.dgf)

# Test corresponding to benchmark C1.2a (constant root params) from Schnepf et al (2020) https://doi.org/10.3389/fpls.2020.00316
# using the numerical method described in Koch et al (2022) https://doi.org/10.1016/j.jcp.2021.110823
dumux_add_test(NAME test_root_soil_benchmark_kernel_c12a
               SOURCES main.cc
               LABELS multidomain multidomain_embedded 1p richards
               TIMEOUT 1500
               CMAKE_GUARD dune-foamgrid_FOUND
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzyData
                        --delimiter " "
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_root_soil_benchmark_kernel_c12a-reference.dat
                                ${CMAKE_CURRENT_BINARY_DIR}/c12a_actualtranspiration.out
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_root_soil_benchmark_kernel_c12a -Root.SpatialParams.Scenario C12a -Output.GnuplotOutputFilenamePrefix c12a_"
)

# Test corresponding to benchmark C1.2b (age-dependent root params) from Schnepf et al (2020) https://doi.org/10.3389/fpls.2020.00316
# using the numerical method described in Koch et al (2022) https://doi.org/10.1016/j.jcp.2021.110823
dumux_add_test(NAME test_root_soil_benchmark_kernel_c12b
               TARGET test_root_soil_benchmark_kernel_c12a
               LABELS multidomain multidomain_embedded 1p richards
               TIMEOUT 1500
               CMAKE_GUARD dune-foamgrid_FOUND
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzyData
                        --delimiter " "
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_root_soil_benchmark_kernel_c12b-reference.dat
                                ${CMAKE_CURRENT_BINARY_DIR}/c12b_actualtranspiration.out
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_root_soil_benchmark_kernel_c12a -Root.SpatialParams.Scenario C12b -Output.GnuplotOutputFilenamePrefix c12b_"
)
