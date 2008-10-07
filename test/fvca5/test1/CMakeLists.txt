add_definitions(-DYASPGRID -DGRIDDIM=2 -DENABLE_UG)

# add build targets
ADD_EXECUTABLE("fvca5_test1" fvca5_test1.cc)
TARGET_LINK_LIBRARIES("fvca5_test1" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
# (this is pretty hacky...)
add_custom_command(TARGET "fvca5_test1"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/../grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/../grids")

