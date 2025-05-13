# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# this is convenience to get a target build_cmake_tests
dune_declare_test_label(LABELS cmake)

# the tests below require dune-common > 2.11
if (dune-common_VERSION VERSION_LESS 2.11)
    return()
endif()

# this to find the dune modules
# in a non-test setup this would be the path to the installed modules
# here we also want to test a build setup without installation
# add the dumux build directory to the path
list(APPEND _DUNE_MODULES_PREFIX_PATH ${PROJECT_BINARY_DIR})
# add the prefix path in case some modules are installed and this is given
list(APPEND _DUNE_MODULES_PREFIX_PATH ${CMAKE_PREFIX_PATH})

# need to escape the semicolon in the path to pass it as a command
string(REPLACE ";" "\;" _DUNE_MODULES_PREFIX_PATH "${_DUNE_MODULES_PREFIX_PATH}")

# this test configures a dumux-depedent-project that doesn't
# use the dune build system at all
add_test(NAME test_dumux_dependent_app_cmake
  COMMAND ${CMAKE_COMMAND}
    -DCMAKE_PREFIX_PATH=${_DUNE_MODULES_PREFIX_PATH}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
    -S ${CMAKE_CURRENT_SOURCE_DIR}/dumux_dependent_project
    -B ${CMAKE_CURRENT_BINARY_DIR}/dumux_dependent_project
)

# this tests building the project
add_test(NAME test_dumux_dependent_app_build
  COMMAND ${CMAKE_COMMAND}
    --build ${CMAKE_CURRENT_BINARY_DIR}/dumux_dependent_project
)

# this tests checks if the executable runs
add_test(NAME test_dumux_dependent_app_test
  COMMAND dumux_dependent_app -NumProcesses 1
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dumux_dependent_project
)

set_tests_properties(test_dumux_dependent_app_build
  PROPERTIES DEPENDS test_dumux_dependent_app_cmake
)

set_tests_properties(test_dumux_dependent_app_test
  PROPERTIES DEPENDS test_dumux_dependent_app_build
)

set_tests_properties(
  test_dumux_dependent_app_build
  test_dumux_dependent_app_cmake
  test_dumux_dependent_app_test
  PROPERTIES LABELS cmake
)

# if the parallel setting verify that we have the same mpi settings
# as dumux. This, for instance, fails if dumux is compiled without MPI,
# but the project has MPI enabled or the other way around.
if (MPI_FOUND)
  add_test(NAME test_dumux_dependent_app_test_parallel
    COMMAND ${MPIEXEC_EXECUTABLE}
      ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 2
      dumux_dependent_app -NumProcesses 2 ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dumux_dependent_project
  )

  set_tests_properties(test_dumux_dependent_app_test_parallel
    PROPERTIES DEPENDS test_dumux_dependent_app_build
  )

  set_tests_properties(test_dumux_dependent_app_test_parallel
    PROPERTIES LABELS cmake
  )
endif()
