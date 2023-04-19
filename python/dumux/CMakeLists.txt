# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(assembly)
add_subdirectory(common)
add_subdirectory(discretization)
add_subdirectory(material)
add_subdirectory(io)
add_subdirectory(porousmediumflow)
add_subdirectory(wrapping)

add_python_targets(dumux
  __init__
)
