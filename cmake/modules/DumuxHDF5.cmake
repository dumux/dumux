# Search HDF5 with C bindings (default)
find_package(HDF5)

if(HDF5_IS_PARALLEL AND HDF5_C_LIBRARIES)
  message(STATUS "Parallel HDF5 C bindings found HDF5_C_LIBRARIES=${HDF5_C_LIBRARIES} HDF5_INCLUDE_DIRS=${HDF5_INCLUDE_DIRS}")
else()
  message(WARNING "Could not find parallel HDF5 C bindings needed for writer. HDF5_IS_PARALLEL=${HDF5_IS_PARALLEL} HDF5_C_LIBRARIES=${HDF5_C_LIBRARIES}"
    " You might need to specify the parallel HDF5 compiler with -DHDF5_C_COMPILER_EXECUTABLE=/usr/bin/h5pcc when calling CMake!")
endif()

if(HDF5_FOUND)
  include_directories(${HDF5_INCLUDE_DIRS})
  link_libraries(${HDF5_C_LIBRARIES})
  add_definitions(${HDF5_DEFINITIONS})
endif(HDF5_FOUND)
