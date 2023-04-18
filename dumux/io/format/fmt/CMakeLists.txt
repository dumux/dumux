# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

if (NOT DUMUX_HAVE_STD_FORMAT)
  dune_library_add_sources(dumux_fmt SOURCES format.cc)

  file(GLOB FMT_HEADERS *.h)
  install(FILES ${FMT_HEADERS} format.cc
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/io/format/fmt)
endif()
