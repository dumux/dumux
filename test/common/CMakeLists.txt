# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(functionfromstringexpression)
add_subdirectory(functions)
add_subdirectory(integrate)
add_subdirectory(math)
add_subdirectory(metadata)
add_subdirectory(numericdifferentiation)
add_subdirectory(parameters)
add_subdirectory(propertysystem)
add_subdirectory(spline)
add_subdirectory(stringutilities)
add_subdirectory(timeloop)
add_subdirectory(typetraits)

dumux_add_test(SOURCES test_partial.cc LABELS unit)
dumux_add_test(SOURCES test_enumerate.cc LABELS unit)
dumux_add_test(SOURCES test_tag.cc LABELS unit)
