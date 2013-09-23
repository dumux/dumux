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
  # create test target for directory, but only if not yet created
  get_directory_test_target(potential_test_target "${CMAKE_CURRENT_BINARY_DIR}")
  if(NOT TARGET ${potential_test_target})
    add_directory_test_target(_test_target)
    set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      PROPERTY TEST_INCLUDE_FILE ${CMAKE_CURRENT_BINARY_DIR}/BuildTests.cmake)
  endif(NOT TARGET ${potential_test_target})

  # add executable
  # check whether executable already exists
  if(NOT TARGET ${dumux_test_executable})
    #set property whether it has to be built with make or only with make test
    if(${DUMUX_BUILD_ALL_TESTS})
      add_executable(${dumux_test_executable} ${dumux_test_executable_source})
    else()
      add_executable(${dumux_test_executable} EXCLUDE_FROM_ALL ${dumux_test_executable_source})
    endif(${DUMUX_BUILD_ALL_TESTS})
  endif(NOT TARGET ${dumux_test_executable})

  # link all libraries to executable
  target_link_libraries(${dumux_test_executable} ${DUNE_LIBS})

  # get optional arguments
  # cannot use ARGN directly with list() command, copy to a variable first
  set(dumux_test_args ${ARGN})
  list(LENGTH dumux_test_args num_dumux_test_args)

  # add test
  add_test(${dumux_test} ${dumux_test_args})
  add_dependencies(${_test_target} ${dumux_test})
endmacro(add_dumux_test)
