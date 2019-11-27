add_input_file_links()

dumux_add_test(NAME test_1pncni_transientbc_box_caching
              LABELS porousmediumflow 1pnc  1pncni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCNITransientBCBox ENABLECACHING=1
              LABELS caching
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2cni_transientbc_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_transientbc_box_caching-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pncni_transientbc_box_caching params.input -Problem.Name test_1p2cni_transientbc_box_caching")

dumux_add_test(NAME test_1pncni_transientbc_box
              LABELS porousmediumflow 1pnc  1pncni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCNITransientBCBox ENABLECACHING=0
              LABELS caching
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2cni_transientbc_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_transientbc_box-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pncni_transientbc_box params.input -Problem.Name test_1p2cni_transientbc_box")

dumux_add_test(NAME test_1pncni_transientbc_tpfa_caching
              LABELS porousmediumflow 1pnc  1pncni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCNITransientBCCCTpfa ENABLECACHING=1
              LABELS caching
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2cni_transientbc_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_transientbc_tpfa_caching-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pncni_transientbc_tpfa_caching params.input -Problem.Name test_1p2cni_transientbc_tpfa_caching")

dumux_add_test(NAME test_1pncni_transientbc_tpfa
              LABELS porousmediumflow 1pnc  1pncni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCNITransientBCCCTpfa ENABLECACHING=0
              LABELS caching
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2cni_transientbc_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_transientbc_tpfa-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pncni_transientbc_tpfa params.input -Problem.Name test_1p2cni_transientbc_tpfa")

dumux_add_test(NAME test_1pncni_transientbc_mpfa_caching
              LABELS porousmediumflow 1pnc  1pncni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCNITransientBCCCMpfa ENABLECACHING=1
              LABELS caching
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2cni_transientbc_mpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_transientbc_mpfa_caching-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pncni_transientbc_mpfa_caching params.input -Problem.Name test_1p2cni_transientbc_mpfa_caching")

dumux_add_test(NAME test_1pncni_transientbc_mpfa
              LABELS porousmediumflow 1pnc  1pncni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCNITransientBCCCMpfa ENABLECACHING=0
              LABELS caching
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2cni_transientbc_mpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_transientbc_mpfa-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pncni_transientbc_mpfa params.input -Problem.Name test_1p2cni_transientbc_mpfa")
