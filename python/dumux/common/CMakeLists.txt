# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
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

target_link_libraries(_common PUBLIC Dune::Common Dune::Grid Dumux::Dumux)

if(SKBUILD)
  install(TARGETS _common LIBRARY DESTINATION python/dumux/common)
endif()
