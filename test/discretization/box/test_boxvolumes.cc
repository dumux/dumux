// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test for finite volume element geometry, sub control volume, and sub
          control volume faces
 */
#include <config.h>

#include <iostream>
#include <utility>

#include <dune/common/test/iteratortest.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/alugrid/grid.hh>
#include <dumux/io/grid/gridmanager_alu.hh>

#include <dumux/common/initialize.hh>
#include <dumux/geometry/volume.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>

int main (int argc, char *argv[])
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);

    std::cout << "Checking the SCVs volumes" << std::endl;

    using Grid = Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>;

    constexpr int dim = Grid::dimension;

    using GridGeometry = BoxFVGridGeometry<double, typename Grid::LeafGridView, true>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GlobalPosition = typename SubControlVolume::GlobalPosition;

    // make a grid
    GlobalPosition lower(0.0);
    GlobalPosition upper(1.0);
    std::array<unsigned int, dim> els{{2, 2, 2}};
    std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, els);
    auto leafGridView = grid->leafGridView();

    GridGeometry gridGeometry(leafGridView);

    for (const auto& element : elements(leafGridView))
    {
        auto eIdx = gridGeometry.elementMapper().index(element);
        auto fvGeometry = localView(gridGeometry);
        fvGeometry.bind(element);

        for (auto&& scv : scvs(fvGeometry))
        {
            const auto& geometry = fvGeometry.geometry(scv);

            const auto idx = scv.dofIndex();
            std::cout << "Volumes for scv: " << "eIdx: " << eIdx << " dofIdx: " << idx << std::endl;
            std::cout << "geometry.volume(): " << geometry.volume() << "   ";
            std::cout << "scv.volume(): " << scv.volume() << "   ";
            std::cout << "Quadrature volume: " << Dumux::volume(geometry) << "   ";
        }
    }
}
