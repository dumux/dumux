# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# Creates symbolic links to all gstat files in the source directory
include_guard(GLOBAL)

macro(add_gstat_file_links)
  FILE(GLOB gstat_files gstat*.txt *.gstat)
  foreach(VAR ${gstat_files})
    get_filename_component(file_name ${VAR} NAME)
    dune_symlink_to_source_files(FILES ${file_name})
  endforeach()
endmacro()
