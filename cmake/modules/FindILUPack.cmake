# .. cmake_module::
#
#    Find the IlUPack library
#
#    You may set the following variables to modify the
#    behaviour of this module:
#
#    :ref:`ILUPACK_ROOT`
#       Path list to search for gstat.
#
#    Sets the following variables:
#
#    :code:`ILUPACK_FOUND`
#       True if the IlUPack library was found.
#
# .. cmake_variable:: ILUPACK_ROOT
#
#   You may set this variable to have :ref:`FindILUPack` look
#   for the gstat library in the given path before inspecting
#   system paths.
#

# look for header files, only at positions given by the user
find_path(ILUPACK_INCLUDE_DIR
    NAMES ilupack.h
    PATHS "${ILUPACK_ROOT}/include"
          "${CMAKE_SOURCE_DIR}/../external/ilupack/"
    PATH_SUFFIXES "include/"
    NO_DEFAULT_PATH
)

# get the root dir of ILUPack
get_filename_component(ILUPACK_DIR ${ILUPACK_INCLUDE_DIR} DIRECTORY)


find_library(LIBILUPACK_MC64
    NAMES libilupack_mc64.a
    PATHS ${ILUPACK_DIR}/lib/GNU64
)

find_library(LIBILUPACK_MUMPS
    NAMES libilupack_mumps.a
    PATHS ${ILUPACK_DIR}/lib/GNU64
)

find_library(LIBGFORTRAN
    NAMES libgfortran.so.3
)

find_library(LIBGOMP
    NAMES libgomp.so.1
)

find_library(LIBMETIS
    NAMES libmetis.a
    PATHS ${ILUPACK_DIR}/lib/GNU64
)

find_library(LIBMUMPS
    NAMES libmumps.a
    PATHS ${ILUPACK_DIR}/lib/GNU64
)

find_library(LIBMETISOMP
    NAMES libmetisomp.a
    PATHS ${ILUPACK_DIR}/lib/GNU64
)

find_library(LIBSPARSPAK
    NAMES libsparspak.a
    PATHS ${ILUPACK_DIR}/lib/GNU64
)

find_library(LIBBLASLIKE
    NAMES libblaslike.a
    PATHS ${ILUPACK_DIR}/lib/GNU64
)

# check version specific macros
include(CheckCSourceCompiles)
include(CMakePushCheckState)
cmake_push_check_state()

# we need if clauses here because variable is set variable-NOTFOUND
#
if(ILUPACK_INCLUDE_DIR)
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${ILUPACK_INCLUDE_DIR})
endif(ILUPACK_INCLUDE_DIR)
if(LIBBLASLIKE)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${LIBBLASLIKE})
endif(LIBBLASLIKE)
if(LIBILUPACK_MUMPS)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${LIBILUPACK_MUMPS})
endif(LIBILUPACK_MUMPS)
if(LIBGFORTRAN)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${LIBGFORTRAN})
endif(LIBGFORTRAN)
if(LIBGOMP)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${LIBGOMP})
endif(LIBGOMP)
if(LIBMUMPS)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${LIBMUMPS})
endif(LIBMUMPS)
if(LIBSPARSPAK)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${LIBSPARSPAK})
endif(LIBSPARSPAK)
if(LIBMETISOMP)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${LIBMETISOMP})
endif(LIBMETISOMP)
if(LIBMETIS)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${LIBMETIS})
endif(LIBMETIS)

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "ILUPACK"
  DEFAULT_MSG
  ILUPACK_INCLUDE_DIR
  LIBBLASLIKE
  LIBILUPACK_MUMPS
  LIBGOMP
  LIBMUMPS
  LIBSPARSPAK
  LIBMETISOMP
  LIBMETIS
)

mark_as_advanced(ILUPACK_INCLUDE_DIR
                 LIBBLASLIKE
                 LIBILUPACK_MUMPS
                 LIBGOMP
                 LIBMUMPS
                 LIBSPARSPAK
                 LIBMETISOMP
                 LIBMETIS
)

# if both headers and libraries are found, store results
if(ILUPACK_FOUND)
    set(ILUPACK_INCLUDE_DIRS ${ILUPACK_INCLUDE_DIR})
    set(ILUPACK_LIBRARIES
        ${LIBBLASLIKE}
        ${LIBILUPACK_MUMPS}
        ${LIBGFORTRAN}
        ${LIBGOMP}
        ${LIBMUMPS}
        ${LIBSPARSPAK}
        ${LIBMETISOMP}
        ${LIBMETIS})
    # log result
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determining location of ILUPack succeeded:\n"
    "Include directory: ${ILUPACK_INCLUDE_DIRS}\n"
    "Library directory: ${ILUPACK_LIBRARIES}\n\n")
    set(ILUPACK_DUNE_COMPILE_FLAGS "${ILUPACK_INCLUDE_DIRS}"
    CACHE STRING "Compile flags used by DUNE when compiling ILUPack programs")
    set(ILUPACK_DUNE_LIBRARIES ${ILUPACK_LIBRARIES}
    CACHE STRING "Libraries used by DUNE when linking ILUPack programs")
else(ILUPACK_FOUND)
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determining location of ILUPack failed:\n"
    "Include directory: ${ILUPACK_INCLUDE_DIRS}\n"
    "Library directory: ${ILUPACK_LIBRARIES}\n\n")
endif(ILUPACK_FOUND)

# set macros for config.h
set(HAVE_ILUPACK ${ILUPACK_FOUND})

if(ILUPACK_FOUND)
    dune_register_package_flags(COMPILE_DEFINITIONS "ENABLE_ILUPACK=1"
                                LIBRARIES "${ILUPACK_LIBRARIES}"
                                INCLUDE_DIRS "${ILUPACK_INCLUDE_DIRS}")
endif()
