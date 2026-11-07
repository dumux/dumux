# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input)

foreach(scheme Box PQ2)
  string(TOLOWER ${scheme} name)
  dumux_add_test(NAME test_multistage_periodic_${name}
                 SOURCES main.cc
                 LABELS unit assembly
                 CMAKE_GUARD "( HAVE_UMFPACK AND dune-spgrid_FOUND )"
                 COMPILE_DEFINITIONS TYPETAG=Heat${scheme}
                 COMMAND ${DUMUX_CONVERGENCETEST}
                 CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_multistage_periodic_${name} params.input -Problem.Name test_multistage_periodic_${name}"
                          --levels "0 true" "1 false" "2 false" "3 false"
                          --level-args "-Grid.Refinement {0} -Test.CompareAssemblers {1}"
                          --errors test_multistage_periodic_${name}-errors.csv
                          --check "L2(u) finest 1.8 2.2")
endforeach()
