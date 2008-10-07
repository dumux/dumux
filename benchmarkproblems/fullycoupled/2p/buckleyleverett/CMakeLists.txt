add_definitions(-DYASPGRID -DGRIDDIM=2 -DENABLE_UG)

# add build targets
ADD_EXECUTABLE("buckleyleverett_fc" buckleyleverett_fc.cc)
TARGET_LINK_LIBRARIES("buckleyleverett_fc" ${DumuxLinkLibraries})

ADD_EXECUTABLE("buckleyleverett1D_fc" buckleyleverett1D_fc.cc)
TARGET_LINK_LIBRARIES("buckleyleverett1D_fc" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
LINK_DIRECTORIES(${DumuxLinkDirectories})
INCLUDE_DIRECTORIES(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "buckleyleverett_fc"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")

