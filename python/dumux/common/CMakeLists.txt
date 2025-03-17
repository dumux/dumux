# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_python_targets(common
  __init__
  boundarytypes
  fvproblem
  fvspatialparams
  properties
)

dune_add_pybind11_module(NAME _common)
target_link_libraries(_common PUBLIC Dune::Common Dune::Grid Dumux::Dumux)

# currently the Python bindings only work with serial dumux code
target_compile_definitions(_common PUBLIC DUMUX_MULTITHREADING_BACKEND=Serial)

if(SKBUILD)
  install(TARGETS _common LIBRARY DESTINATION python/dumux/common)
endif()
