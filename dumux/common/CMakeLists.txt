add_subdirectory(properties)
add_subdirectory(typetraits)

file(GLOB DUMUX_COMMON_HEADERS *.hh *.inc)
install(FILES ${DUMUX_COMMON_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/common)
