# add build targets
add_executable("test_stokes2c" test_stokes2c.cc)
target_link_libraries("test_stokes2c" ${DumuxLinkLibraries})
set_target_properties("test_stokes2c"
  PROPERTIES
  COMPILE_FLAGS ${SUPERLU_CPPFLAGS})

# add required libraries and includes to the build flags
link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories}) 