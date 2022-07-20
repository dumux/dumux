
add_input_file_links()
dune_symlink_to_source_files(FILES grids)

dumux_add_test(NAME test_pnm_2pImbibition_1d_regularized
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1 DIM=2 REGULARIZED=1 DRAINAGE=0
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p-00107.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_1d")

dumux_add_test(NAME test_pnm_2pImbibitionNI_1d_regularized
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=0 DIM=2 REGULARIZED=1 DRAINAGE=0
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p-00107.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_1d")

dumux_add_test(NAME test_pnm_2pImbibition_2d_regularized
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1 DIM=2 REGULARIZED=1 DRAINAGE=0
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p-00107.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_2d")

dumux_add_test(NAME test_pnm_2pImbibition_3d_regularized
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1 DIM=3 REGULARIZED=1 DRAINAGE=0
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p-00107.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_3d")

dumux_add_test(NAME test_pnm_2pImbibition_1d_unregularized
                SOURCES main.cc
                LABELS porenetwork
                COMPILE_DEFINITIONS ISOTHERMAL=1 DIM=2 REGULARIZED=0 DRAINAGE=0
                COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
                CMD_ARGS      --script fuzzy
                                --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p-reference.vtp
                                        ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p-00107.vtp
                                --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_1d")

dumux_add_test(NAME test_pnm_2pImbibitionNI_1d_unregularized
                SOURCES main.cc
                LABELS porenetwork
                COMPILE_DEFINITIONS ISOTHERMAL=0 DIM=2 REGULARIZED=0 DRAINAGE=0
                COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
                CMD_ARGS      --script fuzzy
                                --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p-reference.vtp
                                        ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p-00107.vtp
                                --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_1d")

dumux_add_test(NAME test_pnm_2pImbibition_2d_unregularized
                SOURCES main.cc
                LABELS porenetwork
                COMPILE_DEFINITIONS ISOTHERMAL=1 DIM=2 REGULARIZED=0 DRAINAGE=0
                COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
                CMD_ARGS      --script fuzzy
                                --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p-reference.vtp
                                        ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p-00107.vtp
                                --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_2d")

dumux_add_test(NAME test_pnm_2pImbibition_3d_unregularized
                SOURCES main.cc
                LABELS porenetwork
                COMPILE_DEFINITIONS ISOTHERMAL=1 DIM=3 REGULARIZED=0 DRAINAGE=0
                COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
                CMD_ARGS      --script fuzzy
                                --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p-reference.vtp
                                        ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p-00107.vtp
                                --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_3d")

dumux_add_test(NAME test_pnm_2pImbibitionni
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=0 REGULARIZED=1 DRAINAGE=0
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2pni-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pni-00133.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pni params_ni.input -Problem.Name test_pnm_2pni")
