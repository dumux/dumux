dune_symlink_to_source_files(FILES "params.input" "params_md.input")
dumux_add_test(NAME pnm_3d
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               TIMEOUT "1000"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/pnm_3d params.input")

dumux_add_test(NAME pnm_3d_md
               SOURCES main_md.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1
               COMPILE_DEFINITIONS THROATCONSTRAINT=1
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               TIMEOUT "1000"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/pnm_3d_md params_md.input")
