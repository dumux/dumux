add_subdirectory(dumux)

# link properties.hh needed by the Python bindings
# to determine the list of properties
# create copy for Windows and symlink otherwise
if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
  execute_process(COMMAND ${CMAKE_COMMAND} "-E" "copy" "${CMAKE_SOURCE_DIR}/dumux/common/properties.hh" "${CMAKE_CURRENT_BINARY_DIR}/properties.hh")
else()
  execute_process(COMMAND ${CMAKE_COMMAND} "-E" "create_symlink" "${CMAKE_SOURCE_DIR}/dumux/common/properties.hh" "${CMAKE_CURRENT_BINARY_DIR}/properties.hh")
endif()

dune_python_install_package(
  PATH "."
  PACKAGENAME dumux
  DEPENDS _common
  CMAKE_METADATA_FLAGS DUNE_OPTS_FILE
)
