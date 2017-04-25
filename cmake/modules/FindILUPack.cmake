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
    PATHS "${ILUPACK_ROOT}"
          "${CMAKE_SOURCE_DIR}/../external/ilupack/include/"
)

find_library(LIBILUPACK_MC64
    NAMES libilupack_mc64.a
    PATHS ${CMAKE_SOURCE_DIR}/../external/ilupack/lib/GNU64
)

find_library(LIBILUPACK_MUMPS
    NAMES libilupack_mumps.a
    PATHS ${CMAKE_SOURCE_DIR}/../external/ilupack/lib/GNU64
)

find_library(LIBGFORTRAN
    NAMES libgfortran.so.3
)

find_library(LIBGOMP
    NAMES libgomp.so.1
)

find_library(LIBMETIS
    NAMES libmetis.a
    PATHS ${CMAKE_SOURCE_DIR}/../external/ilupack/lib/GNU64
)

find_library(LIBMUMPS
    NAMES libmumps.a
    PATHS ${CMAKE_SOURCE_DIR}/../external/ilupack/lib/GNU64
)

find_library(LIBMETISOMP
    NAMES libmetisomp.a
    PATHS ${CMAKE_SOURCE_DIR}/../external/ilupack/lib/GNU64
)

find_library(LIBSPARSPAK
    NAMES libsparspak.a
    PATHS ${CMAKE_SOURCE_DIR}/../external/ilupack/lib/GNU64
)

find_library(LIBBLASLIKE
    NAMES libblaslike.a
    PATHS ${CMAKE_SOURCE_DIR}/../external/ilupack/lib/GNU64
)


set(ILUPACK_INCLUDE_DIRS ${ILUPACK_INCLUDE_DIR})
set(ILUPACK_ILUPACK_LIBRARIES
    ${LIBBLASLIKE}
    ${LIBILUPACK_MUMPS}
    ${LIBGFORTRAN}
    ${LIBGOMP}
    ${LIBMUMPS}
    ${LIBSPARSPAK}
    ${LIBMETISOMP}
    ${LIBMETIS})

dune_register_package_flags(COMPILE_DEFINITIONS "ENABLE_ILUPACK=1"
                            LIBRARIES "${ILUPACK_ILUPACK_LIBRARIES}"
                            INCLUDE_DIRS "${ILUPACK_INCLUDE_DIRS}")

# include(FindPackageHandleStandardArgs)
# find_package_handle_standard_args(
#   "ILUPACK"
#   ILUPACK_INCLUDE_DIR
#   #TODO: what does belong here?
# )

# set macros for config.h
set(HAVE_ILUPACK ${ILUPACK_FOUND})
