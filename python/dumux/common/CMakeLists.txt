# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_python_targets(common
  __init__
  boundarytypes
  fvproblem
  fvspatialparams
  properties
)

# currently the Python bindings only work with
# serial dumux code
dune_add_pybind11_module(
  NAME _common
  COMPILE_DEFINITIONS DUMUX_MULTITHREADING_BACKEND=Serial
)

set_property(TARGET _common PROPERTY LINK_LIBRARIES dunecommon dunegrid APPEND)

if(SKBUILD)
  install(TARGETS _common LIBRARY DESTINATION python/dumux/common)
endif()
