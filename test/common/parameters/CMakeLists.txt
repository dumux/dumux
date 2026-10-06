# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(SOURCES test_loggingparametertree.cc
              LABELS unit)
dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(SOURCES test_parametersreset.cc
              LABELS unit)
