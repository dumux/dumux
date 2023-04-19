# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(SOURCES test_propertysystem.cc
              LABELS unit)

dumux_add_test(SOURCES test_propertysystem_aliases.cc
              LABELS unit)

dumux_add_test(SOURCES test_propertysystem_example.cc
              LABELS unit)
