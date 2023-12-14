# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
add_subdirectory(data)
dune_symlink_to_source_files(FILES "params.input")
dune_symlink_to_source_files(FILES "data/matrix_0_0.log.mm")
dune_symlink_to_source_files(FILES "data/matrix_0_1.log.mm")
dune_symlink_to_source_files(FILES "data/matrix_1_0.log.mm")
dune_symlink_to_source_files(FILES "data/matrix_1_1.log.mm")
dune_symlink_to_source_files(FILES "data/rhs_0.log.mm")
dune_symlink_to_source_files(FILES "data/rhs_1.log.mm")

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
dumux_add_test(NAME test_linearsolver_te
               SOURCES test_linearsolver_te.cpp
               LABELS linear unit)
