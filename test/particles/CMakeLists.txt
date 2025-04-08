# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(
  SOURCES test_particles_rotating_disc.cc
  LABELS particles
  COMMAND {CMAKE_CURRENT_BINARY_DIR}/test_particles_rotating_disc
)
