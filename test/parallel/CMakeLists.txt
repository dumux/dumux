# this test only applies if an MPI library is installed!
if (MPI_FOUND)
  add_definitions(-DYASPGRID -DGRIDDIM=2 -DENABLE_UG -DENABLE_MPI)
  add_definitions(${MPI_COMPILE_FLAGS})

  # add build targets
  ADD_EXECUTABLE("test_parallel" test_parallel.cc)
  TARGET_LINK_LIBRARIES("test_parallel" ${DumuxLinkLibraries} ${MPI_LIBRARIES})
  
  # add required libraries and includes to the build flags 
  LINK_DIRECTORIES(${DumuxLinkDirectories})
  INCLUDE_DIRECTORIES(${DumuxIncludeDirectories} ${MPI_INCLUDE_PATH})
  
  # make sure the grids are present in the build directory
  add_custom_command(TARGET "test_parallel"
                     POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E
                          copy_directory 
                             "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                             "${CMAKE_CURRENT_BINARY_DIR}/grids")

endif(MPI_FOUND)
