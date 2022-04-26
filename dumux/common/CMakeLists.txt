add_subdirectory(properties)
add_subdirectory(typetraits)

file(GLOB DUMUX_DUMUX_COMMON_HEADERS *.hh *.inc)
install(FILES ${DUMUX_DUMUX_COMMON_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/dumux/common)
