# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES
    plot_timestepmethods_convergence.py
    plot_timestepmethods_stabilityregions.py
    plot_timestepmethods_symplectic.py
    run_benchmarks.py
)

dumux_add_test(SOURCES test_newmarkbeta.cc LABELS unit timestepping experimental)

dumux_add_test(SOURCES test_timestepmethods.cc LABELS unit timestepping experimental)
set_tests_properties(test_timestepmethods PROPERTIES
  FIXTURES_SETUP timestepmethods_convergence_data)

dumux_add_test(SOURCES test_timestepmethods_stabilityregions.cc LABELS unit timestepping experimental)
set_tests_properties(test_timestepmethods_stabilityregions PROPERTIES
  FIXTURES_SETUP timestepmethods_stability_data)

dumux_add_test(SOURCES test_timestepmethods_symplectic.cc LABELS unit timestepping experimental)
set_tests_properties(test_timestepmethods_symplectic PROPERTIES
  FIXTURES_SETUP timestepmethods_symplectic_data)

dumux_add_test(
    NAME test_timestepmethods_convergence_plot
    TARGET test_timestepmethods
    LABELS unit timestepping experimental
    COMMAND ./plot_timestepmethods_convergence.py
    CMD_ARGS --input test_timestepmethods_data.json
             --output test_timestepmethods_convergence.png
)
set_tests_properties(test_timestepmethods_convergence_plot PROPERTIES
  FIXTURES_REQUIRED timestepmethods_convergence_data)

dumux_add_test(
    NAME test_timestepmethods_stabilityregions_plot
    TARGET test_timestepmethods_stabilityregions
    LABELS unit timestepping experimental
    COMMAND ./plot_timestepmethods_stabilityregions.py
    CMD_ARGS --input test_timestepmethods_stabilityregions_data.json
             --output test_timestepmethods_stabilityregions.png
)
set_tests_properties(test_timestepmethods_stabilityregions_plot PROPERTIES
  FIXTURES_REQUIRED timestepmethods_stability_data)

dumux_add_test(
    NAME test_timestepmethods_symplectic_plot
    TARGET test_timestepmethods_symplectic
    LABELS unit timestepping experimental
    COMMAND ./plot_timestepmethods_symplectic.py
    CMD_ARGS --input test_timestepmethods_symplectic_data.json
             --output test_timestepmethods_symplectic.png
)
set_tests_properties(test_timestepmethods_symplectic_plot PROPERTIES
  FIXTURES_REQUIRED timestepmethods_symplectic_data)
