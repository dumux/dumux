add_subdirectory(doc)
add_subdirectory(img)
add_subdirectory(reference_data)

file(GLOB DUMUX_EXAMPLES_LIDDRIVENCAVITY_HEADERS *.hh *.inc)
install(FILES ${DUMUX_EXAMPLES_LIDDRIVENCAVITY_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/examples/liddrivencavity)
