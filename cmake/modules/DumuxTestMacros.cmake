###
# Add a test. All necessary calls to dune CMake macros for adding a
# test are called from this macro.
# The added test is automatically build if cmake is invode with the command
# line argument -DDUMUX_BUILD_ALL_TESTS:BOOL=TRUE otherwise the test is
# built only when "ctest" / "make test" is invoked.
# The test is only built if it does not already exist.
#
# Arguments:
# - dumux_test:                   name of the new test
# - dumux_test_executable:        name of the executable required by the test
# - dumux_test_executable_source: source file (.cc) of the new test
# - further arguments:            are optional and are used as arguments for calling the test
###
macro(add_dumux_test dumux_test dumux_test_executable dumux_test_executable_source)

  message(WARNING "add_dumux_test is deprecated. Use dune_add_test directly now that we require dune 2.6")

  # if present, symlink the grids folder
  set(grids_directory ${CMAKE_CURRENT_SOURCE_DIR}/grids)
  if(EXISTS ${grids_directory} AND IS_DIRECTORY ${grids_directory})
    dune_symlink_to_source_files(FILES grids)
  endif()

  # if present, symlink the input file
  set(input_file ${CMAKE_CURRENT_SOURCE_DIR}/${dumux_test_executable}.input)
  if(NOT EXISTS ${input_file})
    get_filename_component(base_source_name ${dumux_test_executable_source} NAME_WE)
    set(input_file ${CMAKE_CURRENT_SOURCE_DIR}/${base_source_name}.input)
  endif()
  if(EXISTS ${input_file})
    # dune symlink takes a path relative to the source directory
    get_filename_component(input_file ${input_file} NAME_WE)
    dune_symlink_to_source_files(FILES ${input_file}.input)
  endif()

  # get optional arguments
  # cannot use ARGN directly with list() command, copy to a variable first
  set(dumux_test_args ${ARGN})
  list(LENGTH dumux_test_args num_dumux_test_args)

  # add executable
  # check whether executable already exists
  if(NOT TARGET ${dumux_test_executable})
    add_executable(${dumux_test_executable} ${dumux_test_executable_source})
  endif()

  # add test
  list(GET dumux_test_args 0 dumux_test_command)
  list(REMOVE_AT dumux_test_args 0)
  dune_add_test(NAME ${dumux_test}
                TARGET ${dumux_test_executable}
                COMMAND ${dumux_test_command}
                CMD_ARGS ${dumux_test_args})
  # tests always require the executable to run
  set_tests_properties(${dumux_test} PROPERTIES REQUIRED_FILES ${dumux_test_executable})
endmacro()
