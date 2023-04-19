# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")
dumux_add_test(NAME test_linearsolver
               SOURCES test_linearsolver.cc
               LABELS linear unit)
dumux_add_test(NAME test_parallel_amg_smoothers
               SOURCES test_parallel_amg_smoothers.cc
               LABELS linear unit)
dumux_add_test(NAME test_scalarproducts
               SOURCES test_scalarproducts.cc
               LABELS linear unit)
dumux_add_test(NAME test_dunevectors
               SOURCES test_dunevectors.cc
               LABELS linear unit)
