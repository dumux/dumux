# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input elastica.py run_benchmark.py)

dumux_add_test(NAME test_contact_ring
               LABELS solidmechanics
               SOURCES main.cc
               CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
               COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ringtest.py
               CMD_ARGS test_contact_ring params.input)
