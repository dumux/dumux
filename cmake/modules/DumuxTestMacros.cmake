# Dumux wrapper for the module that provides tools for testing the Dune way.
# We have a wrapper to have to possibily of supporting multiple Dune versions.
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

# Note: This forwards to dune_add_test but enables another layer in case we need to support
# future Dune features with older Dune versions supported by Dumux
function(dumux_add_test)
  dune_add_test(${ARGV})
endfunction()

# Evaluate test guards like dune_add_test internally does
function(dumux_evaluate_cmake_guard GUARD_LETS_YOU_PASS)
  include(CMakeParseArguments)
  set(MULTIARGS CMAKE_GUARD)
  cmake_parse_arguments(EVALGUARD "${OPTIONS}" "${SINGLEARGS}" "${MULTIARGS}" ${ARGN})

  # Check whether the parser produced any errors
  if(EVALGUARD_UNPARSED_ARGUMENTS)
    message(WARNING "Unrecognized arguments ('${EVALGUARD_UNPARSED_ARGUMENTS}') for dumux_evaluate_cmake_guard!")
  endif()

  # determine if all condition of the guard are met
  set(${GUARD_LETS_YOU_PASS} TRUE PARENT_SCOPE)
  set(FAILED_CONDITION_PRINTING "")
  foreach(condition ${EVALGUARD_CMAKE_GUARD})
    separate_arguments(condition)
    if(NOT (${condition}))
      set(${GUARD_LETS_YOU_PASS} FALSE PARENT_SCOPE)
    endif()
  endforeach()
endfunction()
