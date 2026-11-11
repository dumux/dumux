# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

# executable for tpfa tests
dumux_add_test_target(NAME test_1p_convergence_tpfa
                      SOURCES main.cc
                      COMPILE_DEFINITIONS "TYPETAG=OnePIncompressibleTpfa")

# executable for box tests
dumux_add_test_target(NAME test_1p_convergence_box
                      SOURCES main.cc
                      COMPILE_DEFINITIONS "TYPETAG=OnePIncompressibleBox")

# using tpfa and conforming refinement
dumux_add_test(NAME test_1p_convergence_tpfa_conforming
              TARGET test_1p_convergence_tpfa
              CMAKE_GUARD "( dune-functions_FOUND )"
              LABELS porousmediumflow 1p
              COMMAND ./test_1p_convergence_tpfa
              CMD_ARGS params.input -Problem.Name test_1p_convergence_tpfa_conforming)

# using box and conforming refinement
dumux_add_test(NAME test_1p_convergence_box_conforming
              TARGET test_1p_convergence_box
              CMAKE_GUARD "( dune-functions_FOUND )"
              LABELS porousmediumflow 1p
              COMMAND ./test_1p_convergence_box
              CMD_ARGS params.input -Problem.Name test_1p_convergence_box_conforming)
