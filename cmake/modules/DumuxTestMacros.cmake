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
  # if present, copy grids folder
  set(grids_directory ${CMAKE_CURRENT_SOURCE_DIR}/grids)
  if(EXISTS ${grids_directory} AND IS_DIRECTORY ${grids_directory})
    dune_symlink_to_source_files(FILES grids)
  endif()

  # if present, copy input file
  set(input_file ${CMAKE_CURRENT_SOURCE_DIR}/${dumux_test_executable}.input)
  if(NOT EXISTS ${input_file})
    get_filename_component(base_source_name ${dumux_test_executable_source} NAME_WE)
    set(input_file ${CMAKE_CURRENT_SOURCE_DIR}/${base_source_name}.input)
  endif()
  if(EXISTS ${input_file})
    dune_symlink_to_source_files(FILES "${dumux_test_executable}.input")
  endif()

  # add executable
  # check whether executable already exists
  if(NOT TARGET ${dumux_test_executable})
    add_executable(${dumux_test_executable} ${dumux_test_executable_source})
  endif(NOT TARGET ${dumux_test_executable})

  # link all libraries to executable, add all flags
  add_dumux_all_flags(${dumux_test_executable})

  # get optional arguments
  # cannot use ARGN directly with list() command, copy to a variable first
  set(dumux_test_args ${ARGN})
  list(LENGTH dumux_test_args num_dumux_test_args)

  # add test
  add_test(${dumux_test} ${dumux_test_args})

  # return code 77 should be interpreted as skipped test
  set_tests_properties(${dumux_test} PROPERTIES SKIP_RETURN_CODE 77)
endmacro()

###
# adds flags for all third party libraries used in DuMuX to the given targets
###
function(add_dumux_all_flags)
  cmake_parse_arguments(ADD_DUMUX_ALL_FLAGS "" "" "" ${ARGN})
  foreach(_target ${ADD_DUMUX_ALL_FLAGS_UNPARSED_ARGUMENTS})
    # add flags
    add_dune_mpi_flags(${_target})
    add_dune_alugrid_flags(${_target})
    add_dune_parmetis_flags(${_target})
    add_dune_ug_flags(${_target})
    add_dune_umfpack_flags(${_target})
    add_dune_superlu_flags(${_target})
    # add Dune libraries
    target_link_libraries(${_target} ${DUNE_LIBS})
    target_link_libraries(${_target} ${ZLIB_LIBRARIES})
  endforeach()
endfunction()
