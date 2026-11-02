# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# Creates symbolic links to all input files in the source directory
include_guard(GLOBAL)

macro(add_input_file_links)
  message(DEPRECATION "add_input_file_links is deprecated and will be removed after release 3.11. "
                      "Use dune_symlink_to_source_files(FILES ...) and list the input files.")
  FILE(GLOB input_files *.input)
  foreach(VAR ${input_files})
    get_filename_component(file_name ${VAR} NAME)
    dune_symlink_to_source_files(FILES ${file_name})
  endforeach()
endmacro()
