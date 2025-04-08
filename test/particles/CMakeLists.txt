# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(
  SOURCES test_particles_rotating_disc.cc
  LABELS particle
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_particles_rotating_disc
)

dune_symlink_to_source_files(FILES test_particles_fokker_planck.input)
add_executable(test_particles_fokker_planck EXCLUDE_FROM_ALL test_particles_fokker_planck.cc)
foreach(_dim 1 2)
  dumux_add_test(NAME test_particles_fokker_planck_${_dim}d
    TARGET test_particles_fokker_planck
    LABELS particle
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_particles_fokker_planck
    CMD_ARGS -Grid.Dimension ${_dim}
  )
endforeach()
