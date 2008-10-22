add_definitions(-DYASPGRID -DGRIDDIM=2 -DENABLE_UG)

# build target for the simple twophase lens problem
ADD_EXECUTABLE("test_twophase" test_twophase.cc)
TARGET_LINK_LIBRARIES("test_twophase" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "test_twophase"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")

