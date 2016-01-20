# add_dumux_doxgen_target
#
# make sure, that the doxygen links to todo list, bibliography, etc. are correct
MACRO (add_dumux_doxygen_target)
  if(DOXYGEN_FOUND)
    add_doxygen_target()
    add_custom_target(doxygen_${ProjectName}_prebuild
                      COMMAND rm -rf ${CMAKE_BINARY_DIR}/doc/doxygen/html)
    add_dependencies(doxygen_${ProjectName} doxygen_${ProjectName}_prebuild)
    add_custom_command(TARGET doxygen_${ProjectName}
                       POST_BUILD
                       COMMAND ${CMAKE_SOURCE_DIR}/doc/doxygen/sanitizelinks.sh)
  endif()
ENDMACRO ()
