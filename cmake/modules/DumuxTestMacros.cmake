# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# Dumux wrapper for the module that provides tools for testing the Dune way.
# We have a wrapper to have to possibly of supporting multiple Dune versions.
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
#       The numbers of MPI processes that this test should be executed with.
#       Note that one test (in the ctest sense) is created for each number
#       given here. Any number exceeding the user-specified processor maximum
#       :ref:`DUNE_MAX_TEST_CORES` will be ignored. Tests with a
#       processor number :code:`n` higher than one will have the suffix
#       :code:`-mpi-n` appended to their name. You need to specify the
#       TIMEOUT option when specifying the MPI_RANKS option.
#       If the CMD_ARGS contain :code:`--command` (the test runs through
#       :code:`dumux_runtest.py`), the MPI launcher is put in front of the
#       command given after :code:`--command`. Start MPI processes through
#       this option only, so that ctest knows how many processors a test uses.
#
#    .. cmake_param:: THREADS
#       :single:
#
#       The number of threads per process of this test, defaults to
#       :ref:`DUMUX_TEST_NUM_THREADS`. The test runs with the environment variables
#       :code:`DUMUX_NUM_THREADS`, :code:`OMP_NUM_THREADS`, :code:`OPENBLAS_NUM_THREADS`,
#       :code:`MKL_NUM_THREADS` and :code:`VECLIB_MAXIMUM_THREADS` set to this number,
#       and ctest reserves the number of processes times this number of processors for it.
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
#       wrapped around the given commands, unless the CMD_ARGS contain :code:`--command`.
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
#       exists only sometimes.  If your workflow relies on the existence of
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
#    You may set this variable to give an upperbound to the number of MPI processes, that
#    a single test may use. Defaults to 2, when MPI is found and to 1 otherwise.
#    For the tests of DuMux itself, the default is 8.
#
# .. cmake_variable:: DUMUX_TEST_NUM_THREADS
#
#    The number of threads per process of tests that do not specify :code:`THREADS`.
#    Defaults to 1. Set it, e.g. in your opts file, to run all tests multithreaded.
#
# .. cmake_variable:: DUNE_BUILD_TESTS_ON_MAKE_ALL
#
#    You may set this variable through your opts file or on a per module level (in the toplevel
#    :code:`CMakeLists.txt` before :code:`include(DuneMacros)`) to have the Dune build system
#    build all tests during `make all`. Note, that this may take quite some time for some modules.
#    If not in use, you have to build tests through the target :code:`build_tests`.
#
# .. cmake_function:: dumux_evaluate_cmake_guard
#
#    .. cmake_brief::
#
#       Fills the passed variable with TRUE if all guards evaluate to TRUE and FALSE otherwise
#
#    .. cmake_param:: CMAKE_GUARD
#       :multi:
#       :argname: condition
#
#       A number of conditions that CMake should evaluate.
#       Uses the same mechanics that `dumux_add_test` uses to evaluate its CMAKE_GUARD argument.
#
#       The passed condition can be a complex expression like
#       `( A OR B ) AND ( C OR D )`. Mind the spaces around the parentheses.
#
#       Example: Write CMAKE_GUARD dune-foo_FOUND if you want to set a variable
#       that is only true if the module dune-foo has been found.
#
include_guard(GLOBAL)

set(DUMUX_TEST_NUM_THREADS 1 CACHE STRING "Number of threads per process of tests that do not specify THREADS")

# The arguments of the calling dumux_add_test for dune_add_test, without the keywords given to this function
function(_dumux_add_test_forward_args result)
  set(args ${ADDTEST_UNPARSED_ARGUMENTS})
  foreach(keyword IN LISTS OPTIONS)
    if(ADDTEST_${keyword} AND NOT keyword IN_LIST ARGN)
      list(APPEND args ${keyword})
    endif()
  endforeach()
  foreach(keyword IN LISTS SINGLEARGS MULTIARGS)
    if(DEFINED ADDTEST_${keyword} AND NOT keyword IN_LIST ARGN)
      list(APPEND args ${keyword} ${ADDTEST_${keyword}})
    endif()
  endforeach()
  set(${result} ${args} PARENT_SCOPE)
endfunction()

# Note: This forwards to dune_add_test but enables another layer in case we need to support
# future Dune features with older Dune versions supported by Dumux
function(dumux_add_test)
  include(CMakeParseArguments)
  set(OPTIONS EXPECT_COMPILE_FAIL EXPECT_FAIL SKIP_ON_77 COMPILE_ONLY PYTHON_TEST)
  set(SINGLEARGS NAME TARGET TIMEOUT WORKING_DIRECTORY)
  set(MULTIARGS SOURCES COMPILE_DEFINITIONS COMPILE_FLAGS LINK_LIBRARIES CMD_ARGS MPI_RANKS COMMAND CMAKE_GUARD LABELS)
  cmake_parse_arguments(ADDTEST "${OPTIONS}" "${SINGLEARGS};THREADS" "${MULTIARGS}" ${ARGN})

  if(NOT ADDTEST_NAME)
    # try deducing the test name from the executable name
    if(ADDTEST_TARGET)
      set(ADDTEST_NAME ${ADDTEST_TARGET})
    endif()
    # try deducing the test name form the source name
    if(ADDTEST_SOURCES)
      # deducing a name is only possible with a single source argument
      list(LENGTH ADDTEST_SOURCES len)
      if(NOT len STREQUAL "1")
        message(FATAL_ERROR "Cannot deduce test name from multiple sources!")
      endif()
      # strip file extension
      get_filename_component(ADDTEST_NAME ${ADDTEST_SOURCES} NAME_WE)
    endif()
  endif()

  if(NOT DEFINED ADDTEST_THREADS)
    set(ADDTEST_THREADS ${DUMUX_TEST_NUM_THREADS})
  endif()
  if(NOT ADDTEST_THREADS MATCHES "^[1-9][0-9]*$")
    message(FATAL_ERROR "dumux_add_test: THREADS of ${ADDTEST_NAME} has to be a positive integer, not \"${ADDTEST_THREADS}\".")
  endif()

  if(MPIEXEC)
    foreach(arg IN LISTS ADDTEST_COMMAND ADDTEST_CMD_ARGS)
      string(FIND "${arg}" "${MPIEXEC}" position)
      if(NOT position EQUAL -1)
        message(AUTHOR_WARNING "dumux_add_test: ${ADDTEST_NAME} starts ${MPIEXEC} itself, "
                               "so ctest counts it as one process. Use MPI_RANKS instead.")
        break()
      endif()
    endforeach()
  endif()

  dumux_evaluate_cmake_guard(should_run_test CMAKE_GUARD ${ADDTEST_CMAKE_GUARD})
  list(FIND ADDTEST_CMD_ARGS "--command" command_index)
  if(ADDTEST_MPI_RANKS AND NOT command_index EQUAL -1 AND should_run_test AND NOT ADDTEST_COMPILE_ONLY)
    # dune_add_test would start the test script with MPI, so the launcher goes into the command run by the script
    _dumux_add_test_forward_args(build_args CMD_ARGS COMMAND)
    dune_add_test(${build_args} COMPILE_ONLY)
    set(executable ${ADDTEST_TARGET})
    if(ADDTEST_SOURCES)
      set(executable ${ADDTEST_NAME})
    endif()

    set(max_ranks ${DUNE_MAX_TEST_CORES})
    if(NOT HAVE_MPI)
      set(max_ranks 1)
    endif()

    math(EXPR command_index "${command_index} + 1")
    list(GET ADDTEST_CMD_ARGS ${command_index} wrapped_command)
    string(STRIP "${wrapped_command}" wrapped_command)
    if(MPIEXEC_POSTFLAGS)
      list(JOIN MPIEXEC_POSTFLAGS " " postflags)
      string(REGEX REPLACE "^([^ \t\r\n]+)" "\\1 ${postflags}" wrapped_command "${wrapped_command}")
    endif()

    _dumux_add_test_forward_args(test_args NAME TARGET SOURCES COMPILE_DEFINITIONS COMPILE_FLAGS LINK_LIBRARIES CMD_ARGS MPI_RANKS COMPILE_ONLY)
    foreach(procnum ${ADDTEST_MPI_RANKS})
      if(NOT "${procnum}" GREATER "${max_ranks}")
        set(test_name ${ADDTEST_NAME})
        set(cmd_args ${ADDTEST_CMD_ARGS})
        if(NOT ${procnum} STREQUAL "1")
          set(test_name "${test_name}-mpi-${procnum}")
          set(launcher ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${procnum})
          list(JOIN launcher " " launcher)
          list(REMOVE_AT cmd_args ${command_index})
          list(INSERT cmd_args ${command_index} "${launcher} ${wrapped_command}")
        endif()
        dune_add_test(${test_args} NAME ${test_name} TARGET ${executable} CMD_ARGS ${cmd_args})
      endif()
    endforeach()
  else()
    _dumux_add_test_forward_args(forward_args)
    dune_add_test(${forward_args})
  endif()

  if(ADDTEST_SOURCES OR NOT should_run_test)
    set(ADDTEST_TARGET ${ADDTEST_NAME})
  endif()
  target_link_libraries(${ADDTEST_TARGET} PRIVATE Dumux::Dumux)
  dune_target_enable_all_packages(${ADDTEST_TARGET})

  set(thread_environment)
  foreach(variable DUMUX_NUM_THREADS OMP_NUM_THREADS OPENBLAS_NUM_THREADS MKL_NUM_THREADS VECLIB_MAXIMUM_THREADS)
    list(APPEND thread_environment "${variable}=${ADDTEST_THREADS}")
  endforeach()

  if(NOT ADDTEST_MPI_RANKS)
    set(ADDTEST_MPI_RANKS 1)
  endif()

  foreach(procnum ${ADDTEST_MPI_RANKS})
    set(ACTUAL_NAME ${ADDTEST_NAME})
    # add suffix
    if(NOT ${procnum} STREQUAL "1")
      set(ACTUAL_NAME "${ACTUAL_NAME}-mpi-${procnum}")
    endif()

    if(TEST ${ACTUAL_NAME})
      set(ranks 1)
      if(should_run_test)
        set(ranks ${procnum})
      endif()
      math(EXPR processors "${ranks} * ${ADDTEST_THREADS}")
      set_tests_properties(${ACTUAL_NAME} PROPERTIES PROCESSORS ${processors})
      set_property(TEST ${ACTUAL_NAME} APPEND PROPERTY ENVIRONMENT ${thread_environment})

      file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/TestMetaData")
      file(WRITE "${CMAKE_BINARY_DIR}/TestMetaData/${ACTUAL_NAME}.json"
             "{\n  \"name\": \"${ACTUAL_NAME}\",\n  \"target\": \"${ADDTEST_TARGET}\",\n  \"source_dir\": \"${CMAKE_CURRENT_SOURCE_DIR}\"\n}\n")
    endif()
  endforeach()
endfunction()

# Evaluate test guards like dune_add_test internally does
function(dumux_evaluate_cmake_guard GUARD_LETS_YOU_PASS)
  include(CMakeParseArguments)
  set(MULTIARGS CMAKE_GUARD)
  cmake_parse_arguments(EVALGUARD "" "" "${MULTIARGS}" ${ARGN})

  # Check whether the parser produced any errors
  if(EVALGUARD_UNPARSED_ARGUMENTS)
    message(WARNING "Unrecognized arguments ('${EVALGUARD_UNPARSED_ARGUMENTS}') for dumux_evaluate_cmake_guard!")
  endif()

  # determine if all condition of the guard are met
  set(${GUARD_LETS_YOU_PASS} TRUE PARENT_SCOPE)
  foreach(condition ${EVALGUARD_CMAKE_GUARD})
    separate_arguments(condition)
    if(NOT (${condition}))
      set(${GUARD_LETS_YOU_PASS} FALSE PARENT_SCOPE)
    endif()
  endforeach()
endfunction()
