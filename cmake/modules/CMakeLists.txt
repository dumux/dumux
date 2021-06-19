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
  FindJsonCpp.cmake
  FindNLOPT.cmake
  FindPVPython.cmake)
include(GNUInstallDirs)
install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
