add_executable("test_box2pni" test_box2pni.cc)
target_link_libraries("test_box2pni" ${DumuxLinkLibraries})

add_executable("test_cc2pni" test_cc2pni.cc)
target_link_libraries("test_cc2pni" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories})


# make sure the grids are present in the build directory
add_custom_command(TARGET "test_box2pni"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")
