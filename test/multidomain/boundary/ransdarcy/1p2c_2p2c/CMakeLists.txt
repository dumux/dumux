add_input_file_links()
dune_symlink_to_source_files(FILES scripts)

dune_add_test(NAME test_md_boundary_darcy2p2cni_rans1p2cnikepsilon
              SOURCES main.cc
              COMPILE_DEFINITIONS RANSTYPETAG=RANSSubProblemKEpsilonNCNI
              LABELS freeflow rans multidomain
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2cni_rans1p2cnikepsilon_rans-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_rans1p2cnikepsilon_rans-00044.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2cni_rans1p2cnikepsilon_darcy-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_rans1p2cnikepsilon_darcy-00044.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_rans1p2cnikepsilon params.input
                        -Vtk.OutputName test_md_boundary_darcy2p2cni_rans1p2cnikepsilon")

dune_add_test(NAME test_md_boundary_darcy2p2cni_rans1p2cnikomega
              SOURCES main.cc
              COMPILE_DEFINITIONS RANSTYPETAG=RANSSubProblemKOmegaNCNI
              LABELS freeflow rans multidomain
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2cni_rans1p2cnikomega_rans-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_rans1p2cnikomega_rans-00044.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2cni_rans1p2cnikomega_darcy-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_rans1p2cnikomega_darcy-00044.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_rans1p2cnikomega params.input
                        -Vtk.OutputName test_md_boundary_darcy2p2cni_rans1p2cnikomega")

dune_add_test(NAME test_md_boundary_darcy2p2cni_rans1p2cnilowrekepsilon
              SOURCES main.cc
              COMPILE_DEFINITIONS RANSTYPETAG=RANSSubProblemLowReKEpsilonNCNI
              LABELS freeflow rans multidomain
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2cni_rans1p2cnilowrekepsilon_rans-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_rans1p2cnilowrekepsilon_rans-00041.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2cni_rans1p2cnilowrekepsilon_darcy-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_rans1p2cnilowrekepsilon_darcy-00041.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_rans1p2cnilowrekepsilon params.input
                        -Vtk.OutputName test_md_boundary_darcy2p2cni_rans1p2cnilowrekepsilon")

dune_add_test(NAME test_md_boundary_darcy2p2cni_rans1p2cnioneeq
              SOURCES main.cc
              COMPILE_DEFINITIONS RANSTYPETAG=RANSSubProblemOneEqNCNI
              LABELS freeflow rans multidomain
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2cni_rans1p2cnioneeq_rans-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_rans1p2cnioneeq_rans-00041.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2cni_rans1p2cnioneeq_darcy-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_rans1p2cnioneeq_darcy-00041.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_rans1p2cnioneeq params.input
                        -Vtk.OutputName test_md_boundary_darcy2p2cni_rans1p2cnioneeq")

dune_add_test(NAME test_md_boundary_darcy2p2cni_rans1p2cnizeroeq
              SOURCES main.cc
              COMPILE_DEFINITIONS RANSTYPETAG=RANSSubProblemZeroEqNCNI
              LABELS freeflow rans multidomain
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2cni_rans1p2cnizeroeq_rans-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_rans1p2cnizeroeq_rans-00040.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2cni_rans1p2cnizeroeq_darcy-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_rans1p2cnizeroeq_darcy-00040.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_rans1p2cnizeroeq params.input
                        -Vtk.OutputName test_md_boundary_darcy2p2cni_rans1p2cnizeroeq")
