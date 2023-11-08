add_input_file_links()
dune_symlink_to_source_files(FILES grids plot_sw_over_time.py)

dumux_add_test(NAME utrecht_pnm_2pnc
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1 REGULARIZATIONWITHPRESSURE=0
               COMPILE_DEFINITIONS DRAINAGE=1
               COMPILE_DEFINITIONS ONEDDRAINGE=0
               COMPILE_DEFINITIONS CIRCLETHROAT=0
               COMPILE_DEFINITIONS NOREGULARIZTAION=1
               COMPILE_DEFINITIONS ALLOWBREAKTHROUGH=1
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p-00107.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pnc")

dumux_add_test(NAME utrecht_pnm_2pnc_reg
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1 REGULARIZATIONWITHPRESSURE=0
               COMPILE_DEFINITIONS DRAINAGE=1
               COMPILE_DEFINITIONS ONEDDRAINGE=0
               COMPILE_DEFINITIONS CIRCLETHROAT=0
               COMPILE_DEFINITIONS NOREGULARIZTAION=0
               COMPILE_DEFINITIONS ALLOWBREAKTHROUGH=0
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p-00107.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pnc")

dumux_add_test(NAME utrecht_pnm_2pnc_evaporation_reg
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1 REGULARIZATIONWITHPRESSURE=0
               COMPILE_DEFINITIONS DRAINAGE=1
               COMPILE_DEFINITIONS ONEDDRAINGE=0
               COMPILE_DEFINITIONS CIRCLETHROAT=0
               COMPILE_DEFINITIONS NOREGULARIZTAION=0
               COMPILE_DEFINITIONS ALLOWBREAKTHROUGH=1
               COMPILE_DEFINITIONS EVAPORATION=1
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p-00107.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pnc")

dumux_add_test(NAME utrecht_pnm_2pnc_evaporation_reg_circle
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1 REGULARIZATIONWITHPRESSURE=0
               COMPILE_DEFINITIONS DRAINAGE=1
               COMPILE_DEFINITIONS ONEDDRAINGE=0
               COMPILE_DEFINITIONS CIRCLETHROAT=1
               COMPILE_DEFINITIONS NOREGULARIZTAION=0
               COMPILE_DEFINITIONS ALLOWBREAKTHROUGH=0
               COMPILE_DEFINITIONS EVAPORATION=1
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p-00107.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pnc")

dumux_add_test(NAME utrecht_pnm_2pnc_evaporation_circle
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1 REGULARIZATIONWITHPRESSURE=0
               COMPILE_DEFINITIONS DRAINAGE=1
               COMPILE_DEFINITIONS ONEDDRAINGE=0
               COMPILE_DEFINITIONS CIRCLETHROAT=1
               COMPILE_DEFINITIONS NOREGULARIZTAION=1
               COMPILE_DEFINITIONS ALLOWBREAKTHROUGH=0
               COMPILE_DEFINITIONS EVAPORATION=1
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p-00107.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pnc")

dumux_add_test(NAME utrecht_pnm_2pnc_ni
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=0 REGULARIZATIONWITHPRESSURE=0
               COMPILE_DEFINITIONS DRAINAGE=1
               COMPILE_DEFINITIONS ONEDDRAINGE=0
               COMPILE_DEFINITIONS CIRCLETHROAT=0
               COMPILE_DEFINITIONS NOREGULARIZTAION=1
               COMPILE_DEFINITIONS ALLOWBREAKTHROUGH=1
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2pni-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pni-00133.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pnc_ni params_ni.input -Problem.Name test_pnm_2pnc_ni")


dumux_add_test(NAME utrecht_pnm_2pnc_ni_reg
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=0 REGULARIZATIONWITHPRESSURE=0
               COMPILE_DEFINITIONS DRAINAGE=1
               COMPILE_DEFINITIONS ONEDDRAINGE=0
               COMPILE_DEFINITIONS CIRCLETHROAT=0
               COMPILE_DEFINITIONS NOREGULARIZTAION=0
               COMPILE_DEFINITIONS ALLOWBREAKTHROUGH=1
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2pni-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pni-00133.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pnc_ni params_ni.input -Problem.Name test_pnm_2pnc_ni")
