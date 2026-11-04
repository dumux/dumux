# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_spatialparams_ip
              LABELS porousmediumflow unit
              SOURCES test_spatialparams_ip.cc
              COMPILE_DEFINITIONS IPFUNCTION=0)

dumux_add_test(NAME test_spatialparams_ip_fluidmatrixinteraction
              LABELS porousmediumflow unit
              SOURCES test_spatialparams_ip.cc
              COMPILE_DEFINITIONS IPFUNCTION=1
              EXPECT_COMPILE_FAIL)

dumux_add_test(NAME test_spatialparams_ip_wettingphase
              LABELS porousmediumflow unit
              SOURCES test_spatialparams_ip.cc
              COMPILE_DEFINITIONS IPFUNCTION=2
              EXPECT_COMPILE_FAIL)

dumux_add_test(NAME test_spatialparams_ip_characteristiclength
              LABELS porousmediumflow unit
              SOURCES test_spatialparams_ip.cc
              COMPILE_DEFINITIONS IPFUNCTION=3
              EXPECT_COMPILE_FAIL)

dumux_add_test(NAME test_spatialparams_ip_factorenergytransfer
              LABELS porousmediumflow unit
              SOURCES test_spatialparams_ip.cc
              COMPILE_DEFINITIONS IPFUNCTION=4
              EXPECT_COMPILE_FAIL)

dumux_add_test(NAME test_spatialparams_ip_factormasstransfer
              LABELS porousmediumflow unit
              SOURCES test_spatialparams_ip.cc
              COMPILE_DEFINITIONS IPFUNCTION=5
              EXPECT_COMPILE_FAIL)

dumux_add_test(NAME test_spatialparams_ip_inertvolumefraction
              LABELS porousmediumflow unit
              SOURCES test_spatialparams_ip.cc
              COMPILE_DEFINITIONS IPFUNCTION=6
              EXPECT_COMPILE_FAIL)
