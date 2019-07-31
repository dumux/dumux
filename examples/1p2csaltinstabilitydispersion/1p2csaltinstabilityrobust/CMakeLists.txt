add_input_file_links()
dune_symlink_to_source_files(FILES)

#set (CMAKE_BUILD_TYPE Debug)

dune_add_test(NAME 1p2csaltinstabilityrobust
              SOURCES 1p2csaltinstabilityrobust.cc
              COMMAND ${dumux_INCLUDE_DIRS}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/lecture/references/Henry1p2c-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/Henry1p2c-00048.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/1p2csaltinstability -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/1p2csaltinstability.input -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/henry.dgf -TimeLoop.TEnd 1e10")

# headers for installation and headercheck
install(FILES
        1p2csaltinstabilityproblem.hh
        1p2csaltinstabilityspatialparams.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lecture/mm/henryproblem/1p2csaltinstability)
