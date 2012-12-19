if (MPI_FOUND)
  # tell UG that the model is parallelized
  add_definitions(-DModelP)
endif (MPI_FOUND)

# build target for the twophase-twocomponent test problem
add_executable("test_boxrichards" test_boxrichards.cc)
target_link_libraries("test_boxrichards" ${DumuxLinkLibraries})

add_executable("test_ccrichards" test_ccrichards.cc)
target_link_libraries("test_ccrichards" ${DumuxLinkLibraries})

# add required libraries and includes to the build flags 
link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "test_boxrichards"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")

