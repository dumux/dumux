file(GLOB DUMUX_COMMON_PROPERTIES_HEADERS *.hh *.inc)
install(FILES ${DUMUX_COMMON_PROPERTIES_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/common/properties)
