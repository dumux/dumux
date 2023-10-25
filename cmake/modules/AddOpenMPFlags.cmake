# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

include_guard(GLOBAL)

# set variable for config.h
set(DUMUX_HAVE_OPENMP ${OpenMP_FOUND})

# perform DUNE-specific setup tasks
if (OpenMP_FOUND)
  dune_register_package_flags(
    COMPILE_DEFINITIONS ENABLE_OPENMP=1
    LIBRARIES OpenMP::OpenMP_CXX
  )
endif()
