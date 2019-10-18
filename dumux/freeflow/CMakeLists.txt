add_subdirectory(compositional)
add_subdirectory(navierstokes)
add_subdirectory(nonisothermal)
add_subdirectory(rans)
add_subdirectory(shallowwater)

install(FILES
properties.hh
staggeredupwindmethods.hh
turbulencemodel.hh
turbulenceproperties.hh
volumevariables.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow)
