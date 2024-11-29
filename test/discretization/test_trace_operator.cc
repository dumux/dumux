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
#include <dumux/io/grid/facetgridmanager.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/discretization/facetgridmapper.hh>
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
    using TraceGrid = Dune::FoamGrid<1, 2>;
    using TraceGridGeometry = CCTpfaFVGridGeometry<typename TraceGrid::LeafGridView>;
    using TraceGridMapper = FVFacetGridMapper<typename TraceGrid::LeafGridView, GridGeometry>;

    Grid grid{{1.0, 1.0}, {10, 10}};
    auto gridGeometry = std::make_shared<GridGeometry>(grid.leafGridView());

    FacetGridManager<Grid, TraceGrid> traceGridManager;
    traceGridManager.init(grid, [] (const auto&, const auto& is) { return is.boundary(); });
    const auto& traceGridView = traceGridManager.grid().leafGridView();
    auto traceGridGeometry = std::make_shared<TraceGridGeometry>(traceGridView);
    auto traceGridMapper = std::make_shared<TraceGridMapper>(traceGridView, gridGeometry);
    FVTraceOperator traceOperator{traceGridGeometry, traceGridMapper, [&] (const auto& v) {
        return traceGridManager.hostGridVertex(v);
    }};

    int exitCode = 0;
    const auto handleError = [&] (std::string_view message) {
        std::cout << message << std::endl;
        exitCode += 1;
    };

    {
        std::cout << "Testing scv assembly" << std::endl;
        const auto x = makeSolutionVector(*gridGeometry);
        const auto traceX = traceOperator.assembleScvVariables([&] (const auto& scv) {
            return Dune::FieldVector<double, 2>({x[scv.dofIndex()], 0.0});
        });
        for (const auto& v : vertices(traceGridView))
            if (Dune::FloatCmp::ne(traceX[traceGridGeometry->vertexMapper().index(v)][0], testField(v.geometry().center())))
                handleError("Unexpected trace value (scv assembly)");
    }

    {
        std::cout << "Testing scvf assembly" << std::endl;
        const auto traceX = traceOperator.assembleScvfVariables([&] (const auto& scvf) {
            return Dune::FieldVector<double, 2>({testField(scvf.center()), 0.0});
        });
        for (const auto& e : elements(traceGridView))
            if (Dune::FloatCmp::ne(traceX[traceGridGeometry->elementMapper().index(e)][0], testField(e.geometry().center())))
                handleError("Unexpected trace value (scvf assembly)");
    }

    return exitCode;
}
