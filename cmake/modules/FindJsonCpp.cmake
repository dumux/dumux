# .. cmake_module::
#
#    Find the JsonCpp library
#
#    You may set the following variables to modify the
#    behaviour of this module:
#
#    :ref:`JSONCPP_INCLUDE_DIR`
#       Path to include dir of library.
#
#    Sets the following variables:
#
#    :code:`JsonCpp_FOUND`
#       True if the JsonCpp library was found.
#
# .. cmake_variable:: JSONCPP_INCLUDE_DIR
#
#   You may set this variable to have :ref:`FinJsonCpp` a look
#   for the JsonCpp library in the given path before inspecting
#   system paths.
#

# Add a feature summary for this package
include(FeatureSummary)
set_package_properties(JsonCpp PROPERTIES
  DESCRIPTION "JSON for Modern C++ library"
  URL "https://github.com/nlohmann/json"
)

# Try to locate the libraries and their headers, using pkg-config hints
find_path(JSONCPP_INCLUDE_DIR "nlohmann/json.hpp"
          HINTS ${CMAKE_SOURCE_DIR} "${CMAKE_SOURCE_DIR}/.."
          PATH_SUFFIXES "include" "json/include/")

# check version of JSONCpp
find_file(JSONCPP_HEADER "nlohmann/json.hpp"
  HINTS ${JSONCPP_INCLUDE_DIR}
  NO_DEFAULT_PATH)
if(JSONCPP_HEADER)
  file(READ "${JSONCPP_HEADER}" jsoncppheader)
  # get version number from defines in header file
  string(REGEX REPLACE ".*#define NLOHMANN_JSON_VERSION_MAJOR[ \t]+([0-9]+).*" "\\1"
    JSONCPP_MAJOR_VERSION  "${jsoncppheader}")
  string(REGEX REPLACE ".*#define NLOHMANN_JSON_VERSION_MINOR[ \t]+([0-9]+).*" "\\1"
    JSONCPP_MINOR_VERSION  "${jsoncppheader}")
  string(REGEX REPLACE ".*#define NLOHMANN_JSON_VERSION_PATCH[ \t]+([0-9]+).*" "\\1"
    JSONCPP_PATCH_VERSION "${jsoncppheader}")
  if(JSONCPP_MAJOR_VERSION GREATER_EQUAL 0)
    set(JSONCPP_VERSION "${JSONCPP_MAJOR_VERSION}")
  endif()
  if (JSONCPP_MINOR_VERSION GREATER_EQUAL 0)
    set(JSONCPP_VERSION "${JSONCPP_VERSION}.${JSONCPP_MINOR_VERSION}")
  endif()
  if (JSONCPP_PATCH_VERSION GREATER_EQUAL 0)
    set(JSONCPP_VERSION "${JSONCPP_VERSION}.${JSONCPP_PATCH_VERSION}")
  endif()
endif ()

# set if components found
if (JSONCPP_INCLUDE_DIR)
  set(JsonCpp_FOUND TRUE)
endif ()

# dependencies between components
if (NOT JsonCpp_FOUND)
  set(JsonCpp_FOUND FALSE)
endif ()

# Remove these variables from cache inspector
mark_as_advanced(JSONCPP_INCLUDE_DIR)

# Behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args("JsonCpp"
  REQUIRED_VARS
    JSONCPP_INCLUDE_DIR
  VERSION_VAR
    JSONCPP_VERSION
)

# Set targets
if(JsonCpp_FOUND AND NOT TARGET nlohmann_json::nlohmann_json)
  add_library(nlohmann_json::nlohmann_json UNKNOWN IMPORTED)
  set_target_properties(nlohmann_json::nlohmann_json PROPERTIES
    IMPORTED_LOCATION ${JSONCPP_INCLUDE_DIR}
    INTERFACE_INCLUDE_DIRECTORIES ${JSONCPP_INCLUDE_DIR}
  )
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determining location of JsonCpp succeeded:\n"
    "Include directory: ${JSONCPP_INCLUDE_DIR}\n")
else()
  # log erroneous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determining location of JsonCpp failed:\n"
    "Include directory: ${JSONCPP_INCLUDE_DIR}\n")
endif()

# register all JasonCpp related flags
if(JsonCpp_FOUND)
  dune_register_package_flags(INCLUDE_DIRS ${JSONCPP_INCLUDE_DIR})
endif()
