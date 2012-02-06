# add build targets
add_executable("test_stokes" test_stokes.cc)
target_link_libraries("test_stokes" ${DumuxLinkLibraries})
set_target_properties("test_stokes"
  PROPERTIES
  COMPILE_FLAGS ${SUPERLU_CPPFLAGS})
target_link_libraries("test_stokes" ${SUPERLU_LIBS})

# add required libraries and includes to the build flags
link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories}) 