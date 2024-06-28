# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links(FILES params.input)

# isothermal tests

add_executable(test_richards_lens_box_printmed EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_richards_lens_box PUBLIC TYPETAG=RichardsBoxDroplet PUBLIC DIFFMETHOD=DiffMethod::numeric)

dumux_add_test(TARGET test_richards_lens_box_printmed
              LABELS porousmediumflow richards
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_richards_lens_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_richards_lens_box-00007.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_richards_lens_box params.input -Problem.Name test_richards_lens_box")
