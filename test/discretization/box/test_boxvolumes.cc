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

    // initialize parameters
    Parameters::init([](auto& p){
        p["Grid.File"] = "hexahedral.dgf";
    });

    std::cout << "Checking the SCVs volumes" << std::endl;

    using Grid = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>;
    using GridGeometry = BoxFVGridGeometry<double, typename Grid::LeafGridView, true>;

    // make a grid
    GridManager<Grid> gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    GridGeometry gridGeometry(leafGridView);

    double maxRelErrorGeo = 0.0;
    double maxRelErrorScv = 0.0;
    double meanRelErrorGeo = 0.0;
    double meanRelErrorScv = 0.0;
    int n = 0;

    for (const auto& element : elements(leafGridView))
    {
        auto eIdx = gridGeometry.elementMapper().index(element);
        auto fvGeometry = localView(gridGeometry);
        fvGeometry.bind(element);

        for (auto&& scv : scvs(fvGeometry))
        {
            const auto& geometry = fvGeometry.geometry(scv);

            const auto idx = scv.dofIndex();
            const auto exactVol = Dumux::volume(geometry);
            const auto relErrorGeo = std::abs(geometry.volume() - exactVol)/exactVol;
            maxRelErrorGeo = std::max(maxRelErrorGeo, relErrorGeo);
            meanRelErrorGeo += relErrorGeo*relErrorGeo;
            const auto relErrorPScv = std::abs(scv.volume() - exactVol)/exactVol;
            maxRelErrorScv = std::max(maxRelErrorScv, relErrorPScv);
            meanRelErrorScv += relErrorPScv*relErrorPScv;
            std::cout << "Volumes for scv: " << "eIdx: " << eIdx << " dofIdx: " << idx << std::endl;
            std::cout << "Quadrature volume: " << exactVol << "   ";
            std::cout << "geometry.volume() rel error: " << relErrorGeo << "   ";
            std::cout << "scv.volume() rel error: " << relErrorPScv << "   ";
            n++;
        }
    }
    std::cout<<std::endl;
    std::cout << "Max rel error geometry.volume(): " << maxRelErrorGeo << std::endl;
    std::cout << "Max rel error scv.volume() : " << maxRelErrorScv << std::endl;
    std::cout << "Mean squared error geometry.volume(): " << meanRelErrorGeo/n << std::endl;
    std::cout << "Mean squared error scv.volume() : " << meanRelErrorScv/n << std::endl;
}
