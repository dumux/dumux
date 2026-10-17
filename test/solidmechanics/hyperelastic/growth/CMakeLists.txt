# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES convergencetest.py run_benchmark.py)

dumux_add_test(SOURCES test_growth.cc
               LABELS unit solidmechanics)

add_subdirectory(bending_strip)
add_subdirectory(layered_shell)
add_subdirectory(wrinkling)
