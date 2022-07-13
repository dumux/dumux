add_input_file_links()

dumux_add_test(NAME test_pnm_2pnc
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p-00107.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pnc")

dumux_add_test(NAME test_pnm_2pnc_ni
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=0
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2pni-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pni-00133.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pnc_ni params_ni.input -Problem.Name test_pnm_2pnc_ni")
