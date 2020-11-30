set(modules
  AddGstatFileLinks.cmake
  AddInputFileLinks.cmake
  DumuxDoxygen.cmake
  DumuxMacros.cmake
  DumuxTestMacros.cmake
  FindGLPK.cmake
  FindGstat.cmake
  FindNLOPT.cmake
  FindPVPython.cmake)
include(GNUInstallDirs)
install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
