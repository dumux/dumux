add_executable("test_2pdfm" test_2pdfm.cc)
target_link_libraries("test_2pdfm" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags
if(UG_FOUND)
  set_target_properties("test_2pdfm"
    PROPERTIES
    COMPILE_FLAGS "${UG_CPPFLAGS}")
  target_link_libraries("test_2pdfm"
    ${DumuxLinkLibraries}
    ${UG_LIBS})
endif(UG_FOUND)

link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "test_2pdfm"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")
