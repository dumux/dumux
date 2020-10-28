dune_library_add_sources(dumux_fmt SOURCES format.cc)

file(GLOB FMT_HEADERS *.h)
install(FILES ${FMT_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/io/format/fmt)
