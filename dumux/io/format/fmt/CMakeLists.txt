# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

if (NOT DUMUX_HAVE_STD_FORMAT)
  target_sources(dumux_fmt PRIVATE format.cc)

  file(GLOB FMT_HEADERS *.h)
  install(FILES ${FMT_HEADERS} format.cc
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/io/format/fmt)
endif()
