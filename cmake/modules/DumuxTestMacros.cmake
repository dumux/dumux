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

  message(WARNING "add_dumux_test is deprecated. Use dumux_add_test now that we require dune 2.6")

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
  dumux_add_test(NAME ${dumux_test}
                TARGET ${dumux_test_executable}
                COMMAND ${dumux_test_command}
                CMD_ARGS ${dumux_test_args})
  # tests always require the executable to run
  set_tests_properties(${dumux_test} PROPERTIES REQUIRED_FILES ${dumux_test_executable})
endmacro()

# Dumux wrapper for the module that provides tools for testing the Dune way.
# We have a wrapper to have to possibily of supporting multiple Dune versions.
#
# .. cmake_function:: dumux_declare_test_label
#
#    .. cmake_brief::
#
#       Declare labels for :ref:`dumux_add_test`.
#
#    .. cmake_param:: LABELS
#       :multi:
#
#       The names of labels to declare.  Label names must be nonempty and
#       consist only of alphanumeric characters plus :code:`-` and :code:`_`
#       to make sure it is easy to construct regular expressions from them for
#       :code:`ctest -L ${label_regex}`.
#
#    Labels need to be declared to ensure that the target
#    :code:`build_${label}_tests` exists.  They will normally be declared
#    on-demand by :ref:`dumux_add_test`.  But sometimes it is useful to be able to
#    run :code:`make build_${label}_tests` whether or not any tests with that
#    label exists in a module.  For these cases :ref:`dune_declare_test_label` can
#    be called explicitly.
#
#    The label :code:`quick` is always predeclared.
#
# .. cmake_function:: dumux_add_test
#
#    .. cmake_brief::
#
#       Adds a test to the Dumux testing suite!
#
#    .. cmake_param:: NAME
#       :single:
#
#       The name of the test that should be added. If an executable
#       is also added (by specifying SOURCES), the executable is also
#       named accordingly. If omitted, the name will be deduced from
#       the (single) sources parameter or from the given target. Note
#       that this requires you to take care, that you only use a target
#       or source file for but one such test.
#
#    .. cmake_param:: SOURCES
#       :multi:
#
#       The source files that this test depends on. These are the
#       sources that will be passed to :ref:`add_executable`.
#
#       You *must* specify either :code:`SOURCES` or :code:`TARGET`.
#
#    .. cmake_param:: TARGET
#       :single:
#
#       An executable target which should be used for the test. Use
#       this option over the :code:`SOURCES` parameter if you want to
#       reuse already added targets.
#
#       You *must* specify either :code:`SOURCES` or :code:`TARGET`.
#
#    .. cmake_param:: COMPILE_DEFINITIONS
#       :multi:
#       :argname: def
#
#       A set of compile definitions to add to the target.
#       Only definitions beyond the application of :ref:`add_dune_all_flags`
#       have to be stated.
#       This is only used, if :code:`dumux_add_test` adds the executable itself.
#
#    .. cmake_param:: COMPILE_FLAGS
#       :multi:
#       :argname: flag
#
#       A set of non-definition compile flags to add to the target.
#       Only flags beyond the application of :ref:`add_dune_all_flags`
#       have to be stated.
#       This is only used, if :code:`dumux_add_test` adds the executable itself.
#
#    .. cmake_param:: LINK_LIBRARIES
#       :multi:
#       :argname: lib
#
#       A list of libraries to link the target to.
#       Only libraries beyond the application of :ref:`add_dune_all_flags`
#       have to be stated.
#       This is only used, if :code:`dumux_add_test` adds the executable itself.
#
#    .. cmake_param:: EXPECT_COMPILE_FAIL
#       :option:
#
#       If given, the test is expected to not compile successfully!
#
#    .. cmake_param:: EXPECT_FAIL
#       :option:
#
#       If given, this test is expected to compile, but fail to run.
#
#    .. cmake_param:: CMD_ARGS
#       :multi:
#       :argname: arg
#
#       Command line arguments that should be passed to this test.
#
#    .. cmake_param:: MPI_RANKS
#       :multi:
#       :argname: ranks
#
#       The numbers of cores that this test should be executed with.
#       Note that one test (in the ctest sense) is created for each number
#       given here. Any number exceeding the user-specified processor maximum
#       :ref:`DUNE_MAX_TEST_CORES` will be ignored. Tests with a
#       processor number :code:`n` higher than one will have the suffix
#       :code:`-mpi-n` appended to their name. You need to specify the
#       TIMEOUT option when specifying the MPI_RANKS option.
#
#    .. cmake_param:: CMAKE_GUARD
#       :multi:
#       :argname: condition
#
#       A number of conditions that CMake should evaluate before adding this
#       test. If one of the conditions fails, the test should be shown
#       as skipped in the test summary. Use this feature instead of guarding
#       the call to :code:`dumux_add_test` with an :code:`if` clause.
#
#       The passed condition can be a complex expression like
#       `( A OR B ) AND ( C OR D )`. Mind the spaces around the parentheses.
#
#       Example: Write CMAKE_GUARD dune-foo_FOUND if you want your test to only
#       build and run when the dune-foo module is present.
#
#    .. cmake_param:: COMMAND
#       :multi:
#       :argname: cmd
#
#       You may specify the COMMAND option to give the exact command line to be
#       executed when running the test. This defaults to the name of the executable
#       added by dumux_add_test for this test. Note that if you specify both CMD_ARGS
#       and COMMAND, the given CMD_ARGS will be put behind your COMMAND. If you use
#       this in combination with the MPI_RANKS parameter, the call to mpi will still be
#       wrapped around the given commands.
#
#    .. cmake_param:: COMPILE_ONLY
#       :option:
#
#       Set if the given test should only be compiled during :code:`make build_tests`,
#       but not run during :code:`make test`. This is useful if you compile the same
#       executable twice, but with different compile flags, where you want to assure that
#       it compiles with both sets of flags, but you already know they will produce the
#       same result.
#
#    .. cmake_param:: TIMEOUT
#       :single:
#
#       If set, the test will time out after the given number of seconds. This supersedes
#       any timeout setting in ctest (see `cmake --help-property TIMEOUT`). If you
#       specify the MPI_RANKS option, you need to specify a TIMEOUT.
#
#    .. cmake_param:: LABELS
#       :multi:
#
#       A list of labels to add to the test.  This has two effects: it sets
#       the LABELS property on the test so :code:`ctest -L ${label_regex}` can
#       be used to run all tests with certain labels.  It also adds any
#       targets created as dependencies to a custom target, so you can build
#       all tests with a particular label by doing :code:`make
#       build_${label}_tests` without having to build all the other tests as
#       well.
#
#       The :code:`build_${label}_tests` targets are created on-demand the
#       first time a test with that label is added.  In some situations it can
#       depend on the values of cmake cache variables whether a test is added,
#       and then it can happen that the :code:`build_${target}_tests` target
#       exists only sometimes.  If your workflow relies on the existance of
#       these targets, even if building them just returns successfully without
#       doing anything, you can ensure they exist by calling
#       :ref:`dune_declare_test_label` unconditionally.  The label
#       :code:`quick` is always predeclared in this way.
#
#       The label names must be non-empty, and must only contain alphanumeric
#       characters other than :code:`-` or :code:`_`.  This restriction is in
#       place to make it easy to construct regular expressions from the label
#       names for :code:`ctest -L ${label_regex}`.
#
#    This function defines the Dune way of adding a test to the testing suite.
#    You may either add the executable yourself through :ref:`add_executable`
#    and pass it to the :code:`TARGET` option, or you may rely on :ref:`dumux_add_test`
#    to do so.
#
# .. cmake_variable:: DUNE_REENABLE_ADD_TEST
#
#    You may set this variable to True either through your opts file or in your module
#    (before the call to :code:`include(DuneMacros)`) to suppress the error that is thrown if
#    :code:`add_test` is used. You should only do that if you have proper reason to do so.
#
# .. cmake_variable:: DUNE_MAX_TEST_CORES
#
#    You may set this variable to give an upperbound to the number of processors, that
#    a single test may use. Defaults to 2, when MPI is found and to 1 otherwise.
#
# .. cmake_variable:: DUNE_BUILD_TESTS_ON_MAKE_ALL
#
#    You may set this variable through your opts file or on a per module level (in the toplevel
#    :code:`CMakeLists.txt` before :code:`include(DuneMacros)`) to have the Dune build system
#    build all tests during `make all`. Note, that this may take quite some time for some modules.
#    If not in use, you have to build tests through the target :code:`build_tests`.
#

# Note: This is a copy of dune_declare_test_label to be backwards compatible with Dune 2.6 but enable labels
function(dumux_declare_test_label)
  include(CMakeParseArguments)
  set(OPTIONS)
  set(SINGLEARGS)
  set(MULTIARGS LABELS)
  cmake_parse_arguments(arg "${OPTIONS}" "${SINGLEARGS}" "${MULTIARGS}" ${ARGN})

  if( (DEFINED arg_UNPARSED_ARGUMENTS) AND NOT ( arg_UNPARSED_ARGUMENTS STREQUAL "" ) )
    message(FATAL_ERROR "Unhandled extra arguments given to dumux_declare_test_label(): "
      "<${arg_UNPARSED_ARGUMENTS}>")
  endif()

  foreach(label IN LISTS arg_LABELS)
    # Make sure the label is not empty, and does not contain any funny
    # characters, in particular regex characters
    if(NOT (label MATCHES "[-_0-9a-zA-Z]+"))
      message(FATAL_ERROR "Refusing to add label \"${label}\" since it is "
        "empty or contains funny characters (characters other than "
        "alphanumeric ones and \"-\" or \"_\"; the intent of this restriction "
        "is to make construction of the argument to \"ctest -L\" easier")
    endif()
    set(target "build_${label}_tests")
    if(NOT TARGET "${target}")
      add_custom_target("${target}")
    endif()
  endforeach()
endfunction(dumux_declare_test_label)

# predefine "quick" test label so build_quick_tests can be built
# unconditionally
dumux_declare_test_label(LABELS quick)

# Note: This is a copy of dune_declare_test_label to be backwards compatible with Dune 2.6 but enable labels
# After labels are available on a release branch this can simply forward to dune_add_test
function(dumux_add_test)
  # for new versions just forward to dune_add_test
  if(DUNE_COMMON_VERSION VERSION_GREATER 2.6.0)
    dune_add_test(${ARGV})

  # otherwise deal with labels separately (backwards-compatibilty layer with Dune 2.6.0)
  else()
    include(CMakeParseArguments)
    set(OPTIONS EXPECT_COMPILE_FAIL EXPECT_FAIL SKIP_ON_77 COMPILE_ONLY)
    set(SINGLEARGS NAME TARGET TIMEOUT)
    set(MULTIARGS SOURCES COMPILE_DEFINITIONS COMPILE_FLAGS LINK_LIBRARIES CMD_ARGS MPI_RANKS COMMAND CMAKE_GUARD LABELS)
    cmake_parse_arguments(ADDTEST "${OPTIONS}" "${SINGLEARGS}" "${MULTIARGS}" ${ARGN})

    # Check whether the parser produced any errors
    if(ADDTEST_UNPARSED_ARGUMENTS)
      message(WARNING "Unrecognized arguments ('${ADDTEST_UNPARSED_ARGUMENTS}') for dumux_add_test!")
    endif()

    # remove labels from the argument list
    set(FORWARD_ARGS ${ARGV})
    if(ADDTEST_LABELS)
      string(REPLACE "LABELS;${ADDTEST_LABELS}" "" FORWARD_ARGS "${FORWARD_ARGS}")
      # replace head or trailing ";"
      string(REGEX REPLACE ";^" "" FORWARD_ARGS "${FORWARD_ARGS}")
      string(REGEX REPLACE ";$" "" FORWARD_ARGS "${FORWARD_ARGS}")
    endif()

    # foward to dune function
    dune_add_test(${FORWARD_ARGS})

    # take care of labels afterwards
    if(NOT ADDTEST_NAME)
      # try deducing the test name from the executable name
      if(ADDTEST_TARGET)
        set(ADDTEST_NAME ${ADDTEST_TARGET})
      endif()
      # try deducing the test name form the source name
      if(ADDTEST_SOURCES)
        list(LENGTH ADDTEST_SOURCES len)
        get_filename_component(ADDTEST_NAME ${ADDTEST_SOURCES} NAME_WE)
      endif()
    endif()

    if(ADDTEST_SOURCES)
      set(ADDTEST_TARGET ${ADDTEST_NAME})
    endif()

    # possibly set default for mpi ranks
    if(NOT ADDTEST_MPI_RANKS)
      set(ADDTEST_MPI_RANKS 1)
    endif()

    # Discard all parallel tests if MPI was not found
    if(NOT MPI_FOUND)
      set(DUNE_MAX_TEST_CORES 1)
    endif()

    # make sure each label exists and its name is acceptable
    dumux_declare_test_label(LABELS ${ADDTEST_LABELS})

    # Have build_${label}_tests depend on the given target in
    # order to trigger the build correctly
    if(NOT ADDTEST_EXPECT_COMPILE_FAIL)
      foreach(label IN LISTS ADDTEST_LABELS)
        add_dependencies(build_${label}_tests ${ADDTEST_TARGET})
      endforeach()
    endif()

    # Add one test for each specified processor number
    foreach(procnum ${ADDTEST_MPI_RANKS})
      if((NOT "${procnum}" GREATER "${DUNE_MAX_TEST_CORES}") AND (NOT ADDTEST_COMPILE_ONLY))
        set(ACTUAL_NAME ${ADDTEST_NAME})

        if(NOT ${procnum} STREQUAL "1")
          set(ACTUAL_NAME "${ACTUAL_NAME}-mpi-${procnum}")
        endif()

        # Set the labels on the test
        set_tests_properties(${ACTUAL_NAME} PROPERTIES LABELS "${ADDTEST_LABELS}")
      endif()
    endforeach()
  endif()
endfunction()
