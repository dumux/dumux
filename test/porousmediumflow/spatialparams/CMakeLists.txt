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

# the compile-fail tests pass only if the build fails with the expected static assertion,
# WILL_FAIL would invert the regular expression match
set_tests_properties(test_spatialparams_ip_fluidmatrixinteraction PROPERTIES
  WILL_FAIL false PASS_REGULAR_EXPRESSION "auto fluidMatrixInteraction\\(const ElementDiscretization& elemDisc")
set_tests_properties(test_spatialparams_ip_wettingphase PROPERTIES
  WILL_FAIL false PASS_REGULAR_EXPRESSION "int wettingPhase\\(const ElementDiscretization& elemDisc")
set_tests_properties(test_spatialparams_ip_characteristiclength PROPERTIES
  WILL_FAIL false PASS_REGULAR_EXPRESSION "Scalar characteristicLength\\(const ElementDiscretization& elemDisc")
set_tests_properties(test_spatialparams_ip_factorenergytransfer PROPERTIES
  WILL_FAIL false PASS_REGULAR_EXPRESSION "Scalar factorEnergyTransfer\\(const ElementDiscretization& elemDisc")
set_tests_properties(test_spatialparams_ip_factormasstransfer PROPERTIES
  WILL_FAIL false PASS_REGULAR_EXPRESSION "Scalar factorMassTransfer\\(const ElementDiscretization& elemDisc")
set_tests_properties(test_spatialparams_ip_inertvolumefraction PROPERTIES
  WILL_FAIL false PASS_REGULAR_EXPRESSION "Scalar inertVolumeFraction\\(const FVElementGeometry& fvGeometry")
