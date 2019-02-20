add_subdirectory(kepsilon)
add_subdirectory(komega)
add_subdirectory(lowrekepsilon)

install(FILES
indices.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/rans/twoeq)
