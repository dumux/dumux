# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

find_package(Eigen3 3.3 NO_MODULE)
find_package(Spectra 1.0 NO_MODULE HINTS "${DUMUX_SPECTRA_DIR}")

dune_symlink_to_source_files(FILES
    params.input
    run_eigenmodes.py
)

dumux_add_test(NAME test_kirchhoff_love_plate_eigenmodes
              LABELS solidmechanics plate eigenmodes
              SOURCES main.cc
              CMAKE_GUARD "( dune-foamgrid_FOUND AND Eigen3_FOUND AND Spectra_FOUND )"
              COMMAND ./test_kirchhoff_love_plate_eigenmodes)

if(Eigen3_FOUND AND Spectra_FOUND)
    target_link_libraries(test_kirchhoff_love_plate_eigenmodes PRIVATE Eigen3::Eigen Spectra::Spectra)
endif()
