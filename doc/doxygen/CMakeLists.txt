# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
add_subdirectory(PUML)
#dune_symlink_to_source_files(FILES "PUML/HLDD.svg")
file(COPY ${CMAKE_SOURCE_DIR}/doc/doxygen/PUML/HLDD.svg DESTINATION ${CMAKE_BINARY_DIR}/doc/doxygen/PUML)
# shortcut for creating the Doxyfile.in and Doxyfile
add_dumux_doxygen_target()
