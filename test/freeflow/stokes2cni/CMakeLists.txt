# add build targets
add_executable("test_stokes2cni" test_stokes2cni.cc)
target_link_libraries("test_stokes2cni" ${DumuxLinkLibraries})
set_target_properties("test_stokes2cni"
  PROPERTIES
  COMPILE_FLAGS ${SUPERLU_CPPFLAGS})
target_link_libraries("test_stokes2cni" ${SUPERLU_LIBS})

# add required libraries and includes to the build flags
link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories}) 