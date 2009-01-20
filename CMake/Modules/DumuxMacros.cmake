#############################################################
# This sets up the DumuxMacros for the current CMake module.
# Call DumuxSweep at the end of your CMake module in order to
# not pullute the namespace with unused variables
#############################################################
macro(DumuxSetup 
        CMakeModuleName
        ModuleName
        Framework)
  #############
  # Set some internal variables which are used within 
  # the current CMake Module
  set(DumuxModule        ${CMakeModuleName})
  set(DumuxModuleName    ${ModuleName})
  set(DumuxFramework     ${Framework})

  set(DumuxLibsFound 1)
  set(DumuxLibraryNames)
  set(DumuxFound     0)

  set(DumuxPathMessage 
"Set the ${DumuxModule}_DIR cmake cache entry to the directory 
where the ${DumuxModuleName} libraries reside. Alternatively you can set
the ${DumuxFramework}_DIR entry where all ${DumuxFramework} sub-modules have been compiled.")

  # Base path to look for libraries and includes
  if(${DumuxModule}_DIR)
    list(APPEND DumuxModulePath ${${DumuxModule}_DIR})
  endif(${DumuxModule}_DIR)
  if(${DumuxFramework}_DIR)
    list(APPEND DumuxModulePath "${${DumuxFramework}_DIR}/${DumuxModuleName}")
  endif(${DumuxFramework}_DIR)

  # Path to look for includes (->DumuxIncludePath) and libraries (-> DumuxLibraryPath)
  foreach(tmp ${DumuxModulePath})
    list(APPEND DumuxIncludePath "${tmp}" "${tmp}/include")
    list(APPEND DumuxLibraryPath "${tmp}" "${tmp}/lib")
  endforeach(tmp)
  list(APPEND DumuxIncludePath "/usr/include" "/usr/local/include")
  list(APPEND DumuxLibraryPath "/usr/lib" "/usr/local/lib")

  set(DumuxLibraries)
  set(DumuxFailedLibraries)
endmacro(DumuxSetup)

#############################################################
# This adds some additional paths to the location where 
# includes and libraries are searched
#############################################################
macro(DumuxAddPathSuffixes 
        IncludeSuffixes
        LibSuffixes)
  foreach(tmp ${DumuxModulePath})
    # deal with the user defined library locations
    foreach(foo ${LibSuffixes})
      list(APPEND DumuxLibraryPath "${tmp}/${foo}")
    endforeach(foo)

    # deal with the user defined include locations
    foreach(foo ${IncludeSuffixes})
      list(APPEND DumuxIncludePath "${tmp}/${foo}")
    endforeach(foo)
  endforeach(tmp)
endmacro(DumuxAddPathSuffixes)
#############################################################
# Find a given library using some reasonable default 
# search paths. Sets Dumux${LibName}_LIBRARY to the location
# where the library was found and extends the DumuxLibraries
# variable.
#############################################################
macro(DumuxFindLibrary LibName)
  set(Lib ${DumuxModule}_${LibName}_LIBRARY)

  find_library(${Lib}
               ${LibName}
               PATHS ${DumuxLibraryPath}
               PATH_SUFFIXES ".libs")

  if(${Lib})
    list(APPEND DumuxLibraries ${${Lib}})
    list(APPEND DumuxLibraryNames ${LibName})
  else(${Lib})
    list(APPEND DumuxFailedLibraries ${LibName})
  endif(${Lib})
endmacro(DumuxFindLibrary)

#############################################################
# Find a given header file using some reasonable default 
# search paths.
#############################################################
macro(DumuxFindExtraIncludeDir VarName HeaderName)
  set(Inc ${DumuxModule}_${VarName}_INCLUDE_DIR)
  find_path(${Inc}
      ${HeaderName}
      PATHS ${DumuxIncludePath})
    if(${Inc})
      list(APPEND ${DumuxModule}_INCLUDE_DIRS ${${Inc}})
      list(APPEND DumuxIncludes ${Inc})
    else(${Inc})
      list(APPEND DumuxFailedIncludes ${HeaderName})
  endif(${Inc})
endmacro(DumuxFindExtraIncludeDir)

macro(DumuxFindIncludeDir HeaderName)
  set(Inc ${DumuxModule}_INCLUDE_DIR)
  find_path(${Inc}
            ${HeaderName}
            PATHS ${DumuxIncludePath})
  if(${Inc})
    list(APPEND ${DumuxModule}_INCLUDE_DIRS "${${Inc}}")
    list(APPEND DumuxIncludes ${Inc})
  else(${Inc})
    list(APPEND DumuxFailedIncludes ${HeaderName})
  endif(${Inc})
endmacro(DumuxFindIncludeDir)

macro(DumuxFindIncludeBaseDir HeaderName DirSuffix)
  set(Inc ${DumuxModule}_INCLUDE_DIR)
  find_path(${Inc}
            ${HeaderName}
            PATHS ${DumuxIncludePath})
  if(${Inc})
    list(APPEND ${DumuxModule}_INCLUDE_DIRS "${${Inc}}/${DirSuffix}")
    list(APPEND DumuxIncludes ${Inc})
  else(${Inc})
    list(APPEND DumuxFailedIncludes ${HeaderName})
  endif(${Inc})
endmacro(DumuxFindIncludeBaseDir)

#############################################################
# Make sure the required libraries were found
#############################################################
macro(DumuxRequiredLibsFound)
  set(DumuxLibsFound 1)
  set(DumuxFailedLibsMessage "Could not find the required libraries ")
  foreach(curLib ${ARGN})
    set(curLibFound 0)
    foreach(tmp ${DumuxLibraryNames})
      if (tmp STREQUAL ${curLib})
        set(curLibFound 1)
      endif (tmp STREQUAL ${curLib})
    endforeach(tmp)
    
    if (NOT curLibFound)
      set(DumuxLibsFound 0)
      set(DumuxFailedLibsMessage "${DumuxFailedLibsMessage} '${curLib}'")
    endif(NOT curLibFound)
  endforeach(curLib)
endmacro(DumuxRequiredLibsFound)

#############################################################
# Make sure the required libraries were found
#############################################################
macro(DumuxIncludeDirsFound)
endmacro(DumuxIncludeDirsFound)

#############################################################
# Make sure everything required was found
#############################################################
macro(DumuxCheckFound)
  # Set the global macros
  set(DumuxFound 0)

  if(DumuxLibsFound AND ${DumuxModule}_INCLUDE_DIR)
    set(DumuxFound 1)
  endif(DumuxLibsFound AND ${DumuxModule}_INCLUDE_DIR)
  set(${DumuxModule}_FOUND ${DumuxFound})
  set(${DumuxModule}_LIBRARIES ${DumuxLibraries})

  # print status message if requested
  if(NOT ${DumuxModule}_FIND_QUIETLY AND DumuxFound)
    message(STATUS "Found ${DumuxModule}")
  endif(NOT ${DumuxModule}_FIND_QUIETLY AND DumuxFound)

  if(NOT DumuxFound AND ${DumuxModule}_FIND_REQUIRED)
    if (DumuxLibsFound)
      message(FATAL_ERROR "${DumuxPathMessage}")
    else (DumuxLibsFound)
      message(FATAL_ERROR "${DumuxPathMessage}
${DumuxFailedLibsMessage}")
    endif( DumuxLibsFound)
  endif(NOT DumuxFound AND ${DumuxModule}_FIND_REQUIRED)
endmacro(DumuxCheckFound)
