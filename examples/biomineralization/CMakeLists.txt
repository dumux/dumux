add_subdirectory(doc)
add_subdirectory(img)
add_subdirectory(material)

file(GLOB DUMUX_EXAMPLES_BIOMINERALIZATION_HEADERS *.hh *.inc)
install(FILES ${DUMUX_EXAMPLES_BIOMINERALIZATION_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/examples/biomineralization)
