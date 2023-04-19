# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(components)
add_subdirectory(fluidsystems)

add_python_targets(material
  __init__
)
