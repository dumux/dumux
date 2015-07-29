# This file is a copy from dune-common to get these features from 2.4 even with Dune 2.3.
#
# This module provides convenience macros to provide files from the source tree in the build tree.
#
# It provides the following macros:
#
#   dune_add_copy_command(filename)
#
# This macro adds a file-copy command.
# The file_name is the name of a file that exists
# in the source tree. This file will be copied
# to the build tree when executing this command.
# Notice that this does not create a top-level
# target. In order to do this you have to additionally
# call add_custom_target(...) with dependency
# on the file.
#
#   dune_add_copy_target(target_name file_name)
#
# This macro adds a file-copy target under given target_name.
# The file_name is the name of a file that exists
# in the source tree. This file will be copied
# to the build tree.
#
#   dune_add_copy_dependency(target file_name)
#
# This macro adds a copy-dependecy to a target
# The file_name is the name of a file that exists
# in the source tree. This file will be copied
# to the build tree.
#
#   dune_symlink_to_source_tree([NAME name])
#
# add a symlink called NAME to all directories in the build tree (defaults to src_dir).
# That symlink points to the corresponding directory in the source tree.
# Call the macro from the toplevel CMakeLists.txt file of your project.
# You can also call it from some other directory, creating only symlinks
# in that directory and all directories below. A warning is issued on
# Windows systems.
#
#   dune_symlink_to_source_files(FILES files)
#
# add symlinks to the build tree, which point to files in the source tree.
# Foreach file given in "files", a symlink of that name is created in the
# corresponding build directory. Use for ini files, grid files etc. A warning
# is issued on Windows systems.

macro(dune_add_copy_command file_name)
    add_custom_command(
        OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${file_name}"
        COMMAND    ${CMAKE_COMMAND}
        ARGS       -E copy "${CMAKE_CURRENT_SOURCE_DIR}/${file_name}" "${file_name}"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${file_name}"
        )
endmacro(dune_add_copy_command file_name)

macro(dune_add_copy_target target_name file_name)
    dune_add_copy_command(${file_name})
    add_custom_target("${target_name}" ALL DEPENDS "${file_name}")
endmacro(dune_add_copy_target target_name file_name)

macro(dune_add_copy_dependency target file_name)
    message(STATUS "Adding copy-to-build-dir dependency for ${file_name} to target ${target}")
    dune_add_copy_target("${target}_copy_${file_name}" "${file_name}")
    add_dependencies(${target} "${target}_copy_${file_name}")
endmacro(dune_add_copy_dependency)

function(dune_symlink_to_source_tree)
  # parse arguments
  include(CMakeParseArguments)
  cmake_parse_arguments(ARG "" "NAME" "" ${ARGN})
  if(NOT ARG_NAME)
    set(ARG_NAME "src_dir")
  endif()

  # check for Windows to issue a warning
  if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    if(NOT DEFINED DUNE_WINDOWS_SYMLINK_WARNING)
      message(WARNING "Your module wanted to create symlinks, but you cannot do that on your platform.")
      set(DUNE_WINDOWS_SYMLINK_WARNING)
    endif()
  else()
    # get a list of all files in the current source directory and below.
    file(GLOB_RECURSE files RELATIVE ${CMAKE_SOURCE_DIR} "*CMakeLists.txt")

    # iterate over all files, extract the directory name and write a symlink in the corresponding build directory
    foreach(f ${files})
      get_filename_component(dir ${f} DIRECTORY)
      if(NOT "${CMAKE_SOURCE_DIR}/${dir}" MATCHES "${CMAKE_BINARY_DIR}/*")
        execute_process(COMMAND ${CMAKE_COMMAND} "-E" "create_symlink" "${CMAKE_SOURCE_DIR}/${dir}" "${CMAKE_BINARY_DIR}/${dir}/${ARG_NAME}")
      endif()
    endforeach()
  endif()
endfunction(dune_symlink_to_source_tree)

function(dune_symlink_to_source_files)
  # parse arguments
  include(CMakeParseArguments)
  cmake_parse_arguments(ARG "" "" "FILES" ${ARGN})
  if(ARG_UNPARSED_ARGUMENTS)
    message(WARNING "You are using dune_symlink_to_source_files without named arguments (or have typos in your named arguments)!")
  endif()

  # create symlinks for all given files
  foreach(f ${ARG_FILES})
    # check for Windows to issue a warning
    if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
      if(NOT DEFINED DUNE_WINDOWS_SYMLINK_WARNING)
        message(WARNING "Your module wanted to create symlinks, but you cannot do that on your platform.")
        set(DUNE_WINDOWS_SYMLINK_WARNING)
      endif()
      dune_add_copy_command(${f})
    else()
      # create symlink
      execute_process(COMMAND ${CMAKE_COMMAND} "-E" "create_symlink" "${CMAKE_CURRENT_SOURCE_DIR}/${f}" "${CMAKE_CURRENT_BINARY_DIR}/${f}")
    endif()
  endforeach()
endfunction(dune_symlink_to_source_files)
