# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES plot_convergence.py regenerate_doc_images.py disk.geo)

# meshes of the convergence tests, characteristic length clmax
if(Gmsh_FOUND)
  foreach(clmax 0.25 0.125 0.0625)
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/disk_${clmax}.msh
                       COMMAND ${GMSH_EXECUTABLE} -2 -format msh2 -clmax ${clmax} ${CMAKE_CURRENT_SOURCE_DIR}/disk.geo -o ${CMAKE_CURRENT_BINARY_DIR}/disk_${clmax}.msh
                       DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/disk.geo)
    list(APPEND PLATE_DISK_MESHES ${CMAKE_CURRENT_BINARY_DIR}/disk_${clmax}.msh)
  endforeach()
  add_custom_target(test_plate_disk_meshes DEPENDS ${PLATE_DISK_MESHES})
endif()

add_subdirectory(kirchhoff_love)
add_subdirectory(membrane)
add_subdirectory(mindlin_reissner)
