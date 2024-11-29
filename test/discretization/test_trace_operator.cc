// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test trace grid extraction from finite volume grid geometries.
 *
 */
#include <config.h>

#include <vector>
#include <iostream>

#include <dune/grid/yaspgrid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dune/common/fvector.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/discretization/facetgrid.hh>
#include <dumux/discretization/traceoperator.hh>


double testField(const auto& pos)
{ return pos[0]*pos[1]; }


template<typename GridGeometry>
std::vector<double> makeSolutionVector(const GridGeometry& gg)
{
    std::vector<double> result(gg.numDofs());
    for (const auto& element : elements(gg.gridView()))
        for (const auto& scv : scvs(localView(gg).bindElement(element)))
            result[scv.dofIndex()] = testField(scv.dofPosition());
    return result;
}


int main(int argc, char** argv)
{
    using namespace Dumux;
    initialize(argc, argv);

    using Grid = Dune::YaspGrid<2>;
    using GridView = typename Grid::LeafGridView;
    using GridGeometry = BoxFVGridGeometry<double, GridView>;
    using FacetGridType = Dune::FoamGrid<1, 2>;
    using BoundaryGrid = FVFacetGrid<FacetGridType, GridGeometry>;

    Grid grid{{1.0, 1.0}, {10, 10}};
    auto gridGeometry = std::make_shared<GridGeometry>(grid.leafGridView());
    auto traceGrid = std::make_shared<BoundaryGrid>(makeFVBoundaryGrid<FacetGridType>(gridGeometry));
    FVTraceOperator traceOperator{traceGrid};

    {
        std::cout << "Testing scv assembly" << std::endl;
        const auto x = makeSolutionVector(*gridGeometry);
        const auto traceX = traceOperator.assembleScvVariables([&] (const auto& scv) {
            return Dune::FieldVector<double, 2>({x[scv.dofIndex()], 0.0});
        });
        for (const auto& v : vertices(traceGrid->gridView()))
            if (Dune::FloatCmp::ne(traceX[traceGrid->gridView().indexSet().index(v)][0], testField(v.geometry().center())))
            {
                std::cout << "Unexpected trace value" << std::endl;
                return 1;
            }
    }

    {
        std::cout << "Testing scvf assembly" << std::endl;
        const auto traceX = traceOperator.assembleScvfVariables([&] (const auto& scvf) {
            return Dune::FieldVector<double, 2>({testField(scvf.center()), 0.0});
        });
        for (const auto& e : elements(traceGrid->gridView()))
            if (Dune::FloatCmp::ne(traceX[traceGrid->gridView().indexSet().index(e)][0], testField(e.geometry().center())))
            {
                std::cout << "Unexpected trace value" << std::endl;
                return 1;
            }
    }

    return exitCode;
}
