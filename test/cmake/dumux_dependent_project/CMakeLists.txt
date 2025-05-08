# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

cmake_minimum_required(VERSION 3.22)
project(dumux-dependent-project LANGUAGES C CXX VERSION 1.0)
find_package(dumux REQUIRED)

add_executable(dumux_dependent_app main.cc)

# we want to test that linking only against Dumux::Dumux is enough
target_link_libraries(dumux_dependent_app PRIVATE Dumux::Dumux)
