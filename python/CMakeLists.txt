add_subdirectory(dumux)
add_subdirectory(dune)

configure_file(setup.py.in setup.py)

# link properties.hh needed by the Python bindings
# to determine the list of properties
# create copy for Windows and symlink otherwise
if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
  execute_process(COMMAND ${CMAKE_COMMAND} "-E" "copy" "${CMAKE_SOURCE_DIR}/dumux/common/properties.hh" "${CMAKE_CURRENT_BINARY_DIR}/properties.hh")
else()
  execute_process(COMMAND ${CMAKE_COMMAND} "-E" "create_symlink" "${CMAKE_SOURCE_DIR}/dumux/common/properties.hh" "${CMAKE_CURRENT_BINARY_DIR}/properties.hh")
endif()

dune_symlink_to_source_files(FILES "../dumux/common/properties.hh" DESTINATION ".")

dune_python_install_package(
  PATH "."
  CMAKE_METADATA_FILE dune/dumux/metadata.cmake
  DEPENDS _common
  CMAKE_METADATA_FLAGS DUNE_OPTS_FILE
)

if(POLICY CMP0087)
  cmake_policy(SET CMP0087 NEW)
endif()
