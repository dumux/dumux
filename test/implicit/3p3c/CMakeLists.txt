add_executable("test_box3p3c" test_box3p3c.cc)
target_link_libraries("test_box3p3c" ${DumuxLinkLibraries})

add_executable("test_cc3p3c" test_cc3p3c.cc)
target_link_libraries("test_cc3p3c" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories})


# make sure the grids are present in the build directory
add_custom_command(TARGET "test_box3p3c"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")
