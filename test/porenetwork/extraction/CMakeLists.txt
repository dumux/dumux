# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(SOURCES test_spherepackinggeometry.cc
               LABELS porenetwork unit)

dumux_add_test(SOURCES test_spherepackingnetwork.cc
               LABELS porenetwork unit
               CMAKE_GUARD "( dune-foamgrid_FOUND AND CGAL_FOUND )")
if(CGAL_FOUND)
  target_link_libraries(test_spherepackingnetwork PRIVATE CGAL::CGAL)
endif()

dumux_add_test(SOURCES test_spherepackingpermeability.cc
               LABELS porenetwork
               CMAKE_GUARD "( dune-foamgrid_FOUND AND CGAL_FOUND AND HAVE_UMFPACK )")
if(CGAL_FOUND)
  target_link_libraries(test_spherepackingpermeability PRIVATE CGAL::CGAL)
endif()
