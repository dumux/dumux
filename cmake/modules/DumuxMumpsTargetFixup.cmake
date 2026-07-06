# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# Repair the imported MUMPS targets.
#
# The MUMPS package config may bake absolute paths to the (open-)MPI build it was
# compiled against. If that MPI has since been upgraded those paths no longer exist
# and break linking for every consumer of the Dumux interface target. This macro
# re-points any non-existent absolute library dependency of the MUMPS targets to the
# current build's MPI::MPI_Fortran (which tracks the active MPI). It is used both when
# configuring Dumux and, via the exported package config, by downstream modules.
macro(dumux_fix_mumps_targets)
  find_package(MPI QUIET COMPONENTS C Fortran)
  foreach(_mumps_target MUMPS::COMMON MUMPS::DMUMPS)
    if(TARGET ${_mumps_target})
      get_target_property(_mumps_deps ${_mumps_target} INTERFACE_LINK_LIBRARIES)
      if(_mumps_deps)
        set(_mumps_fixed "")
        foreach(_dep IN LISTS _mumps_deps)
          if(IS_ABSOLUTE "${_dep}" AND NOT EXISTS "${_dep}")
            list(APPEND _mumps_fixed MPI::MPI_Fortran)
          else()
            list(APPEND _mumps_fixed "${_dep}")
          endif()
        endforeach()
        list(REMOVE_DUPLICATES _mumps_fixed)
        set_target_properties(${_mumps_target} PROPERTIES INTERFACE_LINK_LIBRARIES "${_mumps_fixed}")
      endif()
    endif()
  endforeach()
endmacro()
