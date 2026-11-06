# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input)

dumux_add_test(NAME test_mindlin_reissner_plate
              LABELS solidmechanics plate
              SOURCES main.cc
              CMAKE_GUARD "( dune-foamgrid_FOUND AND Gmsh_FOUND )"
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_mindlin_reissner_plate params.input"
                       --levels 0.25 0.125 0.0625
                       --level-args "-Grid.File ../disk_{}.msh"
                       --errors test_mindlin_reissner_plate-errors.csv
                       --check "L2Rel(w) mean 1.8 2.2")

if(Gmsh_FOUND)
  add_dependencies(test_mindlin_reissner_plate test_plate_disk_meshes)
endif()
