set(modules
  AddGstatFileLinks.cmake
  AddInputFileLinks.cmake
  DumuxDoxygen.cmake
  DumuxMacros.cmake
  DumuxTestMacros.cmake
  FindGstat.cmake
  FindValgrind.cmake)
include(GNUInstallDirs)
install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
