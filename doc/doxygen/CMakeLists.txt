# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
add_subdirectory(puml)

file(COPY ${CMAKE_SOURCE_DIR}/doc/doxygen/puml/hldd.svg DESTINATION ${CMAKE_BINARY_DIR}/doc/doxygen/puml)
# shortcut for creating the Doxyfile.in and Doxyfile
add_dumux_doxygen_target()
