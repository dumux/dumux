add_definitions(-DYASPGRID -DGRIDDIM=2 -DENABLE_UG)

# add build targets
ADD_EXECUTABLE("convectiondiffusion_ff" convectiondiffusion_ff.cc)
TARGET_LINK_LIBRARIES("convectiondiffusion_ff" ${DumuxLinkLibraries})

ADD_EXECUTABLE("convectiondiffusion1D_ff" convectiondiffusion1D_ff.cc)
TARGET_LINK_LIBRARIES("convectiondiffusion1D_ff" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "convectiondiffusion_ff"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")

