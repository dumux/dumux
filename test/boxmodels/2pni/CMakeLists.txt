ADD_EXECUTABLE("test_2pni" test_2pni.cc)
TARGET_LINK_LIBRARIES("test_2pni" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})


# make sure the grids are present in the build directory
add_custom_command(TARGET "test_2pni"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")
