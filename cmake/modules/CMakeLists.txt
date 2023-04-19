# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

set(modules
  AddGstatFileLinks.cmake
  AddInputFileLinks.cmake
  AddKokkosFlags.cmake
  AddOpenMPFlags.cmake
  DumuxDoxygen.cmake
  DumuxMacros.cmake
  DumuxTestMacros.cmake
  FindGLPK.cmake
  FindGmsh.cmake
  FindGstat.cmake
  FindNLOPT.cmake
  FindPVPython.cmake)
include(GNUInstallDirs)
install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
