# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(adaptive)
add_subdirectory(assembly)
add_subdirectory(common)
add_subdirectory(discretization)
add_subdirectory(flux)
add_subdirectory(freeflow)
add_subdirectory(geomechanics)
add_subdirectory(geometry)
add_subdirectory(io)
add_subdirectory(linear)
add_subdirectory(material)
add_subdirectory(multidomain)
add_subdirectory(nonlinear)
add_subdirectory(parallel)
add_subdirectory(porenetwork)
add_subdirectory(porousmediumflow)

# if Python bindings are enabled, include necessary sub directories.
if(DUMUX_ENABLE_PYTHONBINDINGS)
  add_subdirectory(python)
endif()
