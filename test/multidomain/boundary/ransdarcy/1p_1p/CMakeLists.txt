add_input_file_links()

dune_add_test(NAME test_md_boundary_darcy1p_rans1pkepsilon
              SOURCES main.cc
              COMPILE_DEFINITIONS RANSTYPETAG=RANSSubProblemKEpsilon
              CMAKE_GUARD HAVE_UMFPACK
              LABELS freeflow rans multidomain
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_rans1pkepsilon_darcy-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_rans1pkepsilon_darcy-00031.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_rans1pkepsilon_rans-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_rans1pkepsilon_rans-00031.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_rans1pkepsilon params.input
                        -Vtk.OutputName test_md_boundary_darcy1p_rans1pkepsilon")

dune_add_test(NAME test_md_boundary_darcy1p_rans1pkomega
              SOURCES main.cc
              COMPILE_DEFINITIONS RANSTYPETAG=RANSSubProblemKOmega
              LABELS freeflow rans multidomain
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_rans1pkomega_darcy-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_rans1pkomega_darcy-00033.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_rans1pkomega_rans-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_rans1pkomega_rans-00033.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_rans1pkomega params.input
                        -Vtk.OutputName test_md_boundary_darcy1p_rans1pkomega")

dune_add_test(NAME test_md_boundary_darcy1p_rans1plowrekepsilon
              SOURCES main.cc
              COMPILE_DEFINITIONS RANSTYPETAG=RANSSubProblemLowReKEpsilon
              LABELS freeflow rans multidomain
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_rans1plowrekepsilon_darcy-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_rans1plowrekepsilon_darcy-00035.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_rans1plowrekepsilon_rans-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_rans1plowrekepsilon_rans-00035.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_rans1plowrekepsilon params.input
                        -Vtk.OutputName test_md_boundary_darcy1p_rans1plowrekepsilon")

dune_add_test(NAME test_md_boundary_darcy1p_rans1poneeq
              SOURCES main.cc
              COMPILE_DEFINITIONS RANSTYPETAG=RANSSubProblemOneEq
              LABELS freeflow rans multidomain
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_rans1poneeq_darcy-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_rans1poneeq_darcy-00032.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_rans1poneeq_rans-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_rans1poneeq_rans-00032.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_rans1poneeq params.input
                        -Vtk.OutputName test_md_boundary_darcy1p_rans1poneeq")

dune_add_test(NAME test_md_boundary_darcy1p_rans1pzeroeq
              SOURCES main.cc
              COMPILE_DEFINITIONS RANSTYPETAG=RANSSubProblemZeroEq
              LABELS freeflow rans multidomain
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_rans1pzeroeq_darcy-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_rans1pzeroeq_darcy-00029.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_rans1pzeroeq_rans-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_rans1pzeroeq_rans-00029.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_rans1pzeroeq params.input
                        -Vtk.OutputName test_md_boundary_darcy1p_rans1pzeroeq")

dune_add_test(NAME test_md_boundary_forchheimer1p_rans1pzeroeq
              SOURCES main.cc
              COMPILE_DEFINITIONS RANSTYPETAG=RANSSubProblemZeroEq FORCHHEIMER=1
              LABELS freeflow rans multidomain
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_forchheimer1p_rans1pzeroeq_darcy-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_forchheimer1p_rans1pzeroeq_darcy-00029.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_forchheimer1p_rans1pzeroeq_rans-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_forchheimer1p_rans1pzeroeq_rans-00029.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_forchheimer1p_rans1pzeroeq params.input
                        -Vtk.OutputName test_md_boundary_forchheimer1p_rans1pzeroeq")
