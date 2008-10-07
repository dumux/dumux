add_definitions(-DYASPGRID -DGRIDDIM=2 -DENABLE_UG)

# add build targets
ADD_EXECUTABLE("mcwhorter_ff" mcwhorter_ff.cc)
TARGET_LINK_LIBRARIES("mcwhorter_ff" ${DumuxLinkLibraries})

ADD_EXECUTABLE("mcwhorter1D_ff" mcwhorter1D_ff.cc)
TARGET_LINK_LIBRARIES("mcwhorter1D_ff" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "mcwhorter_ff"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")

