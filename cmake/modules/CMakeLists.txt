set(modules
  AddInputFileLinks.cmake
  CheckConstexpr.cmake
  CopyOfDuneSymlinkOrCopy.cmake
  DumuxDoxygen.cmake
  DumuxMacros.cmake
  DumuxTestMacros.cmake)
include(GNUInstallDirs)
install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
