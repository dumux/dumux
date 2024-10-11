# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

include(DuneAddLibrary)

function(dumux_add_library)
  set(DUNE_ADD_LIB_ARGS "${ARGN}")
  # EXPORT_NAME did not exist prior 2.9, so we remove it from arguments
  if(DUNE_COMMON_VERSION VERSION_LESS 2.9)
      cmake_parse_arguments(ARG "" "EXPORT_NAME" "" ${ARGN})
    string(REPLACE "EXPORT_NAME;${ARG_EXPORT_NAME}" "" DUNE_ADD_LIB_ARGS "${DUNE_ADD_LIB_ARGS}")
  endif()

  # since dune 2.10 dune libraries may have a namespace
  if(DUNE_COMMON_VERSION VERSION_GREATER_EQUAL 2.10)
    set(DUNE_ADD_LIB_ARGS "${DUNE_ADD_LIB_ARGS};NAMESPACE DuMux::")
  endif()

  # forward to dune library creation
  dune_add_library(${DUNE_ADD_LIB_ARGS})
endfunction(dumux_add_library)
