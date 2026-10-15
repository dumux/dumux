# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES
    plot_timestepmethods_convergence.py
    plot_timestepmethods_stabilityregions.py
    plot_timestepmethods_symplectic.py
    plot_multistepmethods_convergence.py
    plot_multistepmethods_stability.py
    plot_generalizedalpha.py
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

dumux_add_test(SOURCES test_multistepmethods.cc LABELS unit timestepping experimental)
set_tests_properties(test_multistepmethods PROPERTIES
  FIXTURES_SETUP multistepmethods_convergence_data)

dumux_add_test(SOURCES test_multistepmethods_stability.cc LABELS unit timestepping experimental)
set_tests_properties(test_multistepmethods_stability PROPERTIES
  FIXTURES_SETUP multistepmethods_stability_data)

dumux_add_test(SOURCES test_generalizedalpha.cc LABELS unit timestepping experimental)
set_tests_properties(test_generalizedalpha PROPERTIES
  FIXTURES_SETUP generalizedalpha_data)

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

dumux_add_test(
    NAME test_multistepmethods_convergence_plot
    TARGET test_multistepmethods
    LABELS unit timestepping experimental
    COMMAND ./plot_multistepmethods_convergence.py
    CMD_ARGS --input test_multistepmethods_data.json
             --output test_multistepmethods_convergence.png
)
set_tests_properties(test_multistepmethods_convergence_plot PROPERTIES
  FIXTURES_REQUIRED multistepmethods_convergence_data)

dumux_add_test(
    NAME test_multistepmethods_stability_plot
    TARGET test_multistepmethods_stability
    LABELS unit timestepping experimental
    COMMAND ./plot_multistepmethods_stability.py
    CMD_ARGS --input test_multistepmethods_stability_data.json
             --output test_multistepmethods_stability.png
)
set_tests_properties(test_multistepmethods_stability_plot PROPERTIES
  FIXTURES_REQUIRED multistepmethods_stability_data)

dumux_add_test(
    NAME test_generalizedalpha_plot
    TARGET test_generalizedalpha
    LABELS unit timestepping experimental
    COMMAND ./plot_generalizedalpha.py
    CMD_ARGS --input test_generalizedalpha_data.json
             --output test_generalizedalpha.png
)
set_tests_properties(test_generalizedalpha_plot PROPERTIES
  FIXTURES_REQUIRED generalizedalpha_data)
