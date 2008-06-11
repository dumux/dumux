#############################################################
# This sets up the StupidMacros for the current CMake module.
# Call StupidSweep at the end of your CMake module in order to
# not pullute the namespace with unused variables
#############################################################
macro(StupidSetup 
        CMakeModuleName
        ModuleName
        Framework)
  #############
  # Set some internal variables which are used within 
  # the current CMake Module
  set(StupidModule        ${CMakeModuleName})
  set(StupidModuleName    ${ModuleName})
  set(StupidFramework     ${Framework})

  set(StupidLibsFound 1)
  set(StupidLibraryNames)
  set(StupidFound         0)

  set(StupidPathMessage 
"Set the ${StupidModule}_DIR cmake cache entry to the directory 
where the ${StupidModuleName} libraries reside. Alternatively you can set
the ${StupidFramework}_DIR entry where all ${StupidFramework} sub-modules have been compiled.")

  # Base path to look for libraries and includes
  if(${StupidModule}_DIR)
    list(APPEND StupidModulePath ${${StupidModule}_DIR})
  endif(${StupidModule}_DIR)
  if(${StupidFramework}_DIR)
    list(APPEND StupidModulePath "${${StupidFramework}_DIR}/${StupidModuleName}")
  endif(${StupidFramework}_DIR)

  # Path to look for includes (->StupidIncludePath) and libraries (-> StupidLibraryPath)
  foreach(tmp ${StupidModulePath})
    list(APPEND StupidIncludePath "${tmp}" "${tmp}/include")
    list(APPEND StupidLibraryPath "${tmp}" "${tmp}/lib")
  endforeach(tmp)
  list(APPEND StupidIncludePath "/usr/include" "/usr/local/include")
  list(APPEND StupidLibraryPath "/usr/lib" "/usr/local/lib")

  set(StupidLibraries)
  set(StupidFailedLibraries)
endmacro(StupidSetup)

#############################################################
# This adds some additional paths to the location where 
# includes and libraries are searched
#############################################################
macro(StupidAddPathSuffixes 
        IncludeSuffixes
        LibSuffixes)
  foreach(tmp ${StupidModulePath})
    # deal with the user defined library locations
    foreach(foo ${LibSuffixes})
      list(APPEND StupidLibraryPath "${tmp}/${foo}")
    endforeach(foo)

    # deal with the user defined include locations
    foreach(foo ${IncludeSuffixes})
      list(APPEND StupidIncludePath "${tmp}/${foo}")
    endforeach(foo)
  endforeach(tmp)
endmacro(StupidAddPathSuffixes)
#############################################################
# Find a given library using some reasonable default 
# search paths. Sets Stupid${LibName}_LIBRARY to the location
# where the library was found and extends the StupidLibraries
# variable.
#############################################################
macro(StupidFindLibrary LibName)
  set(Lib ${StupidModule}_${LibName}_LIBRARY)

  find_library(${Lib}
               ${LibName}
               PATHS ${StupidLibraryPath}
               PATH_SUFFIXES ".libs")

  if(${Lib})
    list(APPEND StupidLibraries ${${Lib}})
    list(APPEND StupidLibraryNames ${LibName})
  else(${Lib})
    list(APPEND StupidFailedLibraries ${LibName})
  endif(${Lib})
endmacro(StupidFindLibrary)

#############################################################
# Find a given header file using some reasonable default 
# search paths.
#############################################################
macro(StupidFindExtraIncludeDir VarName HeaderName)
  set(Inc ${StupidModule}_${VarName}_INCLUDE_DIR)
  find_path(${Inc}
      ${HeaderName}
      PATHS ${StupidIncludePath})
    if(${Inc})
      list(APPEND ${StupidModule}_INCLUDE_DIRS ${${Inc}})
      list(APPEND StupidIncludes ${Inc})
    else(${Inc})
      list(APPEND StupidFailedIncludes ${HeaderName})
  endif(${Inc})
endmacro(StupidFindExtraIncludeDir)

macro(StupidFindIncludeDir HeaderName)
  set(Inc ${StupidModule}_INCLUDE_DIR)
  find_path(${Inc}
            ${HeaderName}
            PATHS ${StupidIncludePath})
  if(${Inc})
    list(APPEND ${StupidModule}_INCLUDE_DIRS "${${Inc}}")
    list(APPEND StupidIncludes ${Inc})
  else(${Inc})
    list(APPEND StupidFailedIncludes ${HeaderName})
  endif(${Inc})
endmacro(StupidFindIncludeDir)

macro(StupidFindIncludeBaseDir HeaderName DirSuffix)
  set(Inc ${StupidModule}_INCLUDE_DIR)
  find_path(${Inc}
            ${HeaderName}
            PATHS ${StupidIncludePath})
  if(${Inc})
    list(APPEND ${StupidModule}_INCLUDE_DIRS "${${Inc}}/${DirSuffix}")
    list(APPEND StupidIncludes ${Inc})
  else(${Inc})
    list(APPEND StupidFailedIncludes ${HeaderName})
  endif(${Inc})
endmacro(StupidFindIncludeBaseDir)

#############################################################
# Make sure the required libraries were found
#############################################################
macro(StupidRequiredLibsFound)
  set(StupidLibsFound 1)
  set(StupidFailedLibsMessage "Could not find the required libraries ")
  foreach(curLib ${ARGN})
    set(curLibFound 0)
    foreach(tmp ${StupidLibraryNames})
      if (tmp STREQUAL ${curLib})
        set(curLibFound 1)
      endif (tmp STREQUAL ${curLib})
    endforeach(tmp)
    
    if (NOT curLibFound)
      set(StupidLibsFound 0)
      set(StupidFailedLibsMessage "${StupidFailedLibsMessage} '${curLib}'")
    endif(NOT curLibFound)
  endforeach(curLib)
endmacro(StupidRequiredLibsFound)

#############################################################
# Make sure the required libraries were found
#############################################################
macro(StupidIncludeDirsFound)
endmacro(StupidIncludeDirsFound)

#############################################################
# Make sure everything required was found
#############################################################
macro(StupidCheckFound)
  # Set the global macros
  set(StupidFound 0)

  if(StupidLibsFound AND ${StupidModule}_INCLUDE_DIR)
    set(StupidFound 1)
  endif(StupidLibsFound AND ${StupidModule}_INCLUDE_DIR)
  set(${StupidModule}_FOUND ${StupidFound})
  set(${StupidModule}_LIBRARIES ${StupidLibraries})

  # print status message if requested
  if(NOT ${StupidModule}_FIND_QUIETLY AND StupidFound)
    message(STATUS "Found ${StupidModule}")
  endif(NOT ${StupidModule}_FIND_QUIETLY AND StupidFound)

  if(NOT StupidFound AND ${StupidModule}_FIND_REQUIRED)
    if (StupidLibsFound)
      message(FATAL_ERROR "${StupidPathMessage}")
    else (StupidLibsFound)
      message(FATAL_ERROR "${StupidPathMessage}
${StupidFailedLibsMessage}")
    endif( StupidLibsFound)
  endif(NOT StupidFound AND ${StupidModule}_FIND_REQUIRED)
endmacro(StupidCheckFound)
