# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(SOURCES test_nonlinear_preconditioning_partition.cc
               LABELS unit nonlinear preconditioning)

dumux_add_test(SOURCES test_nonlinear_preconditioning_views.cc
               LABELS unit nonlinear preconditioning)

dumux_add_test(SOURCES test_nonlinear_preconditioning_subdomain_assembly.cc
               LABELS unit nonlinear preconditioning
               CMD_ARGS params.input)

dune_symlink_to_source_files(FILES params.input)

dumux_add_test(SOURCES test_nonlinear_preconditioning_schwarznewton.cc
               LABELS unit nonlinear preconditioning
               CMD_ARGS params.input)

dumux_add_test(SOURCES test_nonlinear_preconditioning_robustness.cc
               LABELS nonlinear preconditioning
               CMD_ARGS params.input -Newton.MaxTimeStepDivisions 1 -Newton.Verbosity 0)

dumux_add_test(SOURCES test_nonlinear_preconditioning_scaling.cc
               LABELS nonlinear preconditioning
               CMD_ARGS params.input -Grid.Cells "32 32" -Newton.Verbosity 0)

dumux_add_test(SOURCES test_nonlinear_preconditioning_pvswitch.cc
               LABELS unit nonlinear preconditioning
               CMD_ARGS params_2p2c.input)

dune_symlink_to_source_files(FILES params_2p2c.input)

dumux_add_test(SOURCES test_nonlinear_preconditioning_richards1d.cc
               LABELS unit nonlinear preconditioning
               CMD_ARGS params_richards1d.input)

dune_symlink_to_source_files(FILES params_richards1d.input)
