add_subdirectory(2pinfiltration)
add_subdirectory(1ptracer)
add_subdirectory(biomineralization)
add_subdirectory(shallowwaterfriction)
add_subdirectory(freeflowchannel)
add_subdirectory(1protationsymmetry)
add_subdirectory(liddrivencavity)

# if Python bindings are enabled, include Python binding tests
if(DUNE_ENABLE_PYTHONBINDINGS)
  add_subdirectory(python)
endif()
