add_definitions(-DYASPGRID -DGRIDDIM=2 -DENABLE_UG)

# add build targets
ADD_EXECUTABLE("test_box3d" test_box3d.cc)
TARGET_LINK_LIBRARIES("test_box3d" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "test_box3d"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")

