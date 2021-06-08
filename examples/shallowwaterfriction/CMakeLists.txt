add_subdirectory(doc)
add_subdirectory(img)

file(GLOB DUMUX_EXAMPLES_SHALLOWWATERFRICTION_HEADERS *.hh *.inc)
install(FILES ${DUMUX_EXAMPLES_SHALLOWWATERFRICTION_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/examples/shallowwaterfriction)
