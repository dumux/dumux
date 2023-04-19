# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input convergencetest.py)

dumux_add_test(NAME test_1p_rootbenchmark_tpfa
               LABELS porousmediumflow 1p
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=RootBenchmarkCCTpfa
               COMMAND ./convergencetest.py
               CMD_ARGS test_1p_rootbenchmark_tpfa -Problem.EnableGravity false)

dumux_add_test(NAME test_1p_rootbenchmark_tpfa_gravity
               TARGET test_1p_rootbenchmark_tpfa
               LABELS porousmediumflow 1p
               COMPILE_DEFINITIONS TYPETAG=RootBenchmarkCCTpfa
               COMMAND ./convergencetest.py
               CMD_ARGS test_1p_rootbenchmark_tpfa test_1p_rootbenchmark_tpfa_gravity.log)
