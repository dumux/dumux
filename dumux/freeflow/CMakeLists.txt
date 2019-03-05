add_subdirectory(compositional)
add_subdirectory(navierstokes)
add_subdirectory(nonisothermal)
add_subdirectory(rans)

install(FILES
properties.hh
turbulenceproperties.hh
volumevariables.hh
staggeredupwindmethods.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow)
