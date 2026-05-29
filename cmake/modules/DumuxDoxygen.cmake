# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# add_dumux_doxgen_target
#
# make sure, that the doxygen links to todo list, bibliography, etc. are correct
include_guard(GLOBAL)

include(DuneDoxygen)

macro (add_dumux_doxygen_target)
  if(DOXYGEN_FOUND)
    add_doxygen_target()
    add_custom_target(doxygen_${ProjectName}_prebuild
                      COMMAND rm -rf ${CMAKE_BINARY_DIR}/doc/doxygen/html)

    if(TARGET doxygen_${ProjectName})
      add_dependencies(doxygen_${ProjectName} doxygen_${ProjectName}_prebuild)
    elseif(TARGET doxygen_${ProjectName}_build)
      add_dependencies(doxygen_${ProjectName}_build doxygen_${ProjectName}_prebuild)
    else()
      message(WARNING
        "Neither doxygen_${ProjectName} nor doxygen_${ProjectName}_build exists; "
        "cannot add DuMux Doxygen prebuild dependency"
      )
    endif()

  endif()
endmacro ()
