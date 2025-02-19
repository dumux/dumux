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
    using TraceGrid = Dune::FoamGrid<1, 2>;
    using TraceGridGeometry = CCTpfaFVGridGeometry<typename TraceGrid::LeafGridView>;

    Grid grid{{1.0, 1.0}, {10, 10}};

    int exitCode = 0;
    const auto handleError = [&] (std::string_view message) {
        std::cout << message << std::endl;
        exitCode += 1;
    };

    {
        std::cout << "Testing box trace operator" << std::endl;
        using GridGeometry = BoxFVGridGeometry<double, GridView>;
        using TraceGridMapper = FVFacetGridMapper<typename TraceGrid::LeafGridView, GridGeometry>;
        auto gridGeometry = std::make_shared<GridGeometry>(grid.leafGridView());

        FacetGridManager<Grid, TraceGrid> traceGridManager;
        traceGridManager.init(grid, [] (const auto&, const auto& is) { return is.boundary(); });
        const auto& traceGridView = traceGridManager.grid().leafGridView();

        auto traceGridGeometry = std::make_shared<TraceGridGeometry>(traceGridView);
        auto traceGridMapper = std::make_shared<TraceGridMapper>(traceGridView, gridGeometry);
        TraceOperatorFactory factory{traceGridGeometry, traceGridMapper, [&] (const auto& v) {
            return traceGridManager.hostGridVertex(v);
        }};

        {
            std::cout << " -- testing trace operator" << std::endl;
            auto traceOperator = factory.traceOperator();
            const auto traceX = traceOperator.apply([&] (const auto& scvf, const auto& context) {
                return Dune::FieldVector<double, 2>({
                    testField(context.gridGeometryLocalView().scv(scvf.insideScvIdx()).dofPosition()),
                    0.0
                });
            });
            if (traceX.size() != traceGridGeometry->vertexMapper().size())
                handleError("Unexpected trace coefficient vector size (box)");
            else
                for (const auto& v : vertices(traceGridView))
                {
                    const auto x = traceX[traceGridGeometry->vertexMapper().index(v)][0];
                    const auto expected = testField(v.geometry().center());
                    if (Dune::FloatCmp::ne(x, expected))
                        handleError("Unexpected trace value (box): " + std::to_string(x) + " vs " + std::to_string(expected));
                }
        }

        {
            std::cout << " -- testing normal trace operator" << std::endl;
            auto traceOperator = factory.normalTraceOperator();
            const auto traceX = traceOperator.apply([&] (const auto& scvf, const auto&) {
                return testField(scvf.center());
            });
            for (const auto& e : elements(traceGridView))
                if (Dune::FloatCmp::ne(traceX[traceGridGeometry->elementMapper().index(e)], testField(e.geometry().center())))
                    handleError("Unexpected normal trace value (box)");
        }
    }

    {
        std::cout << "Testing cc trace operator" << std::endl;
        using GridGeometry = CCTpfaFVGridGeometry<GridView>;
        using TraceGridMapper = FVFacetGridMapper<typename TraceGrid::LeafGridView, GridGeometry>;
        auto gridGeometry = std::make_shared<GridGeometry>(grid.leafGridView());

        FacetGridManager<Grid, TraceGrid> traceGridManager;
        traceGridManager.init(grid, [] (const auto&, const auto& is) { return is.boundary(); });
        const auto& traceGridView = traceGridManager.grid().leafGridView();

        auto traceGridGeometry = std::make_shared<TraceGridGeometry>(traceGridView);
        auto traceGridMapper = std::make_shared<TraceGridMapper>(traceGridView, gridGeometry);
        TraceOperatorFactory factory{traceGridGeometry, traceGridMapper};

        {
            std::cout << " -- testing trace operator" << std::endl;
            auto traceOperator = factory.traceOperator();
            const auto traceX = traceOperator.apply([&] (const auto& scvf, const auto& context) {
                return Dune::FieldVector<double, 2>({
                    testField(scvf.center()),
                    0.0
                });
            });
            if (traceX.size() != traceGridGeometry->elementMapper().size())
                handleError("Unexpected trace coefficient vector size (tpfa)");
            else
                for (const auto& e : elements(traceGridView))
                {
                    const auto x = traceX[traceGridGeometry->elementMapper().index(e)][0];
                    const auto expected = testField(e.geometry().center());
                    if (Dune::FloatCmp::ne(x, expected))
                        handleError("Unexpected trace value (tpfa): " + std::to_string(x) + " vs " + std::to_string(expected));
                }
        }

        {
            std::cout << " -- testing normal trace operator" << std::endl;
            auto traceOperator = factory.normalTraceOperator();
            const auto traceX = traceOperator.apply([&] (const auto& scvf, const auto&) {
                return testField(scvf.center());
            });
            for (const auto& e : elements(traceGridView))
                if (Dune::FloatCmp::ne(traceX[traceGridGeometry->elementMapper().index(e)], testField(e.geometry().center())))
                    handleError("Unexpected normal trace value (tpfa)");
        }
    }

    return exitCode;
}
