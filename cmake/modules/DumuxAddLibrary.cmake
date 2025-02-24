# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# Dumux function to add a library
#
# .. cmake_function:: dumux_add_library
#
#    .. cmake_brief::
#
#       Add a library to the build system. For details see :ref:`dune_add_library`.
#
include(DuneAddLibrary)

function(dumux_add_library)
  set(DUNE_ADD_LIB_ARGS "${ARGN}")

  # EXPORT_NAME and NAMESPACE did not exist prior 2.10, so we remove it from arguments
  if(DUNE_COMMON_VERSION VERSION_LESS 2.10)
    set(SINGLEARGS EXPORT_NAME NAMESPACE)
    cmake_parse_arguments(ARG "" "${SINGLEARGS}" "" ${ARGN})
    string(REPLACE "EXPORT_NAME;${ARG_EXPORT_NAME}" "" DUNE_ADD_LIB_ARGS "${DUNE_ADD_LIB_ARGS}")
    string(REPLACE "NAMESPACE;${ARG_NAMESPACE}" "" DUNE_ADD_LIB_ARGS "${DUNE_ADD_LIB_ARGS}")
  endif()

  # forward to dune library creation
  dune_add_library(${DUNE_ADD_LIB_ARGS})
endfunction()
