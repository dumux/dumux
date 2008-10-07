add_definitions(-DYASPGRID -DGRIDDIM=2 -DENABLE_UG)

# add build targets
ADD_EXECUTABLE("mcwhorter_fc" mcwhorter_fc.cc)
TARGET_LINK_LIBRARIES("mcwhorter_fc" ${DumuxLinkLibraries})

ADD_EXECUTABLE("mcwhorter1D_fc" mcwhorter1D_fc.cc)
TARGET_LINK_LIBRARIES("mcwhorter1D_fc" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "mcwhorter_fc"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")

