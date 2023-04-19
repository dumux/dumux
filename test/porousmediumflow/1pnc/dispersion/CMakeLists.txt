# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

dumux_add_test(NAME test_1pnc_dispersion_scheidegger
              LABELS porousmediumflow dispersion 1pnc
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNCDispersionTestBox SCHEIDEGGER=1
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1pnc_dispersion_scheidegger-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1pnc_dispersion_scheidegger-00024.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pnc_dispersion_scheidegger params.input
                        -Problem.Name test_1pnc_dispersion_scheidegger
                        -Problem.AlphaL 0.1
                        -Problem.AlphaT 0.01")

dumux_add_test(NAME test_1pnc_dispersion_fulldispersiontensor_isotropic
              LABELS porousmediumflow dispersion 1pnc
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNCDispersionTestBox SCHEIDEGGER=0
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1pnc_dispersion_fulldispersiontensor_isotropic-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1pnc_dispersion_fulldispersiontensor_isotropic-00025.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pnc_dispersion_fulldispersiontensor_isotropic params.input
                        -Problem.Name test_1pnc_dispersion_fulldispersiontensor_isotropic
                        -Problem.DispersionTensor 1e-6")

dumux_add_test(NAME test_1pnc_dispersion_fulldispersiontensor_anisotropic
              LABELS porousmediumflow dispersion 1pnc
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNCDispersionTestBox SCHEIDEGGER=0
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1pnc_dispersion_fulldispersiontensor_anisotropic-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1pnc_dispersion_fulldispersiontensor_anisotropic-00025.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pnc_dispersion_fulldispersiontensor_anisotropic params.input
                        -Problem.Name test_1pnc_dispersion_fulldispersiontensor_anisotropic
                        -Problem.DispersionTensor '1e-6 5e-8 6e-9 5e-6'")

dumux_add_test(NAME test_1pncni_dispersion_scheidegger
              LABELS porousmediumflow dispersion 1pncni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNCNIDispersionTestBox SCHEIDEGGER=1 NONISOTHERMAL=1
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1pncni_dispersion_scheidegger-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1pncni_dispersion_scheidegger-00026.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pncni_dispersion_scheidegger params.input
                        -Problem.Name test_1pncni_dispersion_scheidegger
                        -Problem.AlphaL 0.1
                        -Problem.AlphaT 0.01")

dumux_add_test(NAME test_1pncni_dispersion_fulldispersiontensor_isotropic
              LABELS porousmediumflow dispersion 1pncni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNCNIDispersionTestBox SCHEIDEGGER=0 NONISOTHERMAL=1
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1pncni_dispersion_fulldispersiontensor_isotropic-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1pncni_dispersion_fulldispersiontensor_isotropic-00026.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pncni_dispersion_fulldispersiontensor_isotropic params.input
                        -Problem.Name test_1pncni_dispersion_fulldispersiontensor_isotropic
                        -Problem.DispersionTensor 1e-6")

dumux_add_test(NAME test_1pncni_dispersion_fulldispersiontensor_anisotropic
              LABELS porousmediumflow dispersion 1pncni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNCNIDispersionTestBox SCHEIDEGGER=0 NONISOTHERMAL=1
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1pncni_dispersion_fulldispersiontensor_anisotropic-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1pncni_dispersion_fulldispersiontensor_anisotropic-00026.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pncni_dispersion_fulldispersiontensor_anisotropic params.input
                        -Problem.Name test_1pncni_dispersion_fulldispersiontensor_anisotropic
                        -Problem.DispersionTensor '1e-6 5e-8 6e-9 5e-6'")
