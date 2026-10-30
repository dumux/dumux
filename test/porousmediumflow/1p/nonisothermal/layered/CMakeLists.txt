# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input)

dumux_add_test(NAME test_1pni_layered_box
               SOURCES main.cc
               LABELS porousmediumflow 1p 1pni
               COMMAND ./test_1pni_layered_box
               CMD_ARGS params.input)
