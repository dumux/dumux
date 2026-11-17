// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test of the integration point source and its helper
 *
 * The weight of a point source (quadrature weight times integration element)
 * is split among all elements and sub-control volumes that contain the point.
 * Summed over all entries of the point source map, the split weights have to
 * recover the weight of the original point source, also for points that lie on
 * element boundaries and sub-control volume boundaries.
 */
#include <config.h>

#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/method.hh>
#include <dumux/geometry/intersectingentities.hh>
#include <dumux/multidomain/embedded/integrationpointsource.hh>

namespace Dumux {

constexpr int dim = 3;
using GlobalPosition = Dune::FieldVector<double, dim>;
using TestSource = IntegrationPointSource<GlobalPosition, Dune::FieldVector<double, 1>>;

void setLumping(bool lumping)
{
    Parameters::init([&](Dune::ParameterTree& params){
        params["PointSource.EnableBoxLumping"] = lumping ? "true" : "false";
    });
}

template<class GridGeometry>
std::vector<GlobalPosition> testPositions(const GridGeometry& gridGeometry)
{
    const auto& gridView = gridGeometry.gridView();
    std::vector<GlobalPosition> positions;

    // points exactly on the vertices, edges, faces and in the centers of the elements
    for (const auto& element : elements(gridView))
    {
        const auto geometry = element.geometry();
        const auto refElement = Dune::referenceElement(geometry);
        for (int codim = 0; codim <= dim; ++codim)
            for (int i = 0; i < refElement.size(codim); ++i)
                positions.push_back(geometry.global(refElement.position(i, codim)));
    }

    // random points
    std::mt19937 generator(42);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    for (int i = 0; i < 200; ++i)
    {
        GlobalPosition p;
        for (int d = 0; d < dim; ++d)
            p[d] = distribution(generator);
        positions.push_back(p);
    }

    return positions;
}

template<class GridGeometry>
bool testConservation(const GridGeometry& gridGeometry, bool lumping, const std::string& name, std::size_t& numEntries)
{
    std::mt19937 generator(7);
    std::uniform_real_distribution<double> distribution(0.1, 2.0);

    // one point source per element containing the point, as in the coupling managers
    std::vector<TestSource> sources;
    double inputWeight = 0.0;
    std::size_t id = 0;
    const auto& tree = gridGeometry.boundingBoxTree();
    for (const auto& position : testPositions(gridGeometry))
    {
        const auto entities = intersectingEntities(position, tree);
        if (entities.empty())
            continue;

        const auto weight = distribution(generator);
        const auto integrationElement = distribution(generator);
        inputWeight += weight*integrationElement;
        for (const auto eIdx : entities)
        {
            sources.emplace_back(position, id, weight, integrationElement, eIdx);
            sources.back().setEmbeddings(entities.size());
        }
        ++id;
    }

    setLumping(lumping);
    std::map<std::pair<std::size_t, std::size_t>, std::vector<TestSource>> pointSourceMap;
    IntegrationPointSourceHelper::computePointSourceMap(gridGeometry, sources, pointSourceMap);

    double mapWeight = 0.0;
    numEntries = 0;
    for (const auto& [key, entries] : pointSourceMap)
        for (const auto& s : entries)
        {
            ++numEntries;
            mapWeight += s.quadratureWeight()*s.integrationElement()/s.embeddings();
        }

    const auto error = std::abs(mapWeight - inputWeight)/inputWeight;
    std::cout << name << (lumping ? ", lumping: " : ", basis function weights: ")
              << sources.size() << " sources, relative weight error " << error << std::endl;
    return error < 1e-13;
}

template<class Grid>
int testGrid(const Grid& grid, const std::string& gridName)
{
    int failures = 0;
    const auto gridView = grid.leafGridView();

    {
        BoxFVGridGeometry<double, typename Grid::LeafGridView> gridGeometry(gridView);
        std::size_t numEntriesLumping = 0, numEntriesBasisFunction = 0;
        if (!testConservation(gridGeometry, true, gridName + ", box", numEntriesLumping))
            ++failures;
        if (!testConservation(gridGeometry, false, gridName + ", box", numEntriesBasisFunction))
            ++failures;

        // the parameter is read for every call: spreading over all sub-control volumes
        // creates more entries than the assignment to the sub-control volumes containing the point
        if (numEntriesBasisFunction <= numEntriesLumping)
        {
            std::cout << gridName << ": the box lumping parameter has no effect" << std::endl;
            ++failures;
        }
    }

    {
        CCTpfaFVGridGeometry<typename Grid::LeafGridView> gridGeometry(gridView);
        std::size_t numEntries = 0;
        if (!testConservation(gridGeometry, true, gridName + ", cctpfa", numEntries))
            ++failures;
    }

    return failures;
}

bool testSetters()
{
    TestSource source(GlobalPosition(0.5), Dune::FieldVector<double, 1>(1.0), 3, 0.25, 2.0, 5);
    source.setQuadratureWeight(0.75);
    source.setIntegrationElement(4.0);
    source.setEmbeddings(3);
    return source.quadratureWeight() == 0.75
        && source.integrationElement() == 4.0
        && source.embeddings() == 3
        && source.elementIndex() == 5
        && source.id() == 3;
}

// a point source that does not lie in any sub-control volume of its element must not be dropped silently
template<class GridGeometry>
bool testThrowsOutsideElement(const GridGeometry& gridGeometry)
{
    std::map<std::pair<std::size_t, std::size_t>, std::vector<TestSource>> pointSourceMap;
    const std::vector<TestSource> sources{ TestSource(GlobalPosition(10.0), 0, 1.0, 1.0, 0) };
    setLumping(true);
    try { IntegrationPointSourceHelper::computePointSourceMap(gridGeometry, sources, pointSourceMap); }
    catch (const Dune::InvalidStateException&) { return true; }
    return false;
}

struct PQ1BubbleGridGeometry
{ static constexpr auto discMethod = DiscretizationMethods::pq1bubble; };

// discretization methods the point source map is not implemented for must throw
bool testThrowsUnsupportedMethod()
{
    std::map<std::pair<std::size_t, std::size_t>, std::vector<TestSource>> pointSourceMap;
    const std::vector<TestSource> sources{ TestSource(GlobalPosition(0.5), 0, 1.0, 1.0, 0) };
    try { IntegrationPointSourceHelper::computePointSourceMap(PQ1BubbleGridGeometry{}, sources, pointSourceMap); }
    catch (const Dune::NotImplemented&) { return true; }
    return false;
}

} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;
    initialize(argc, argv);
    Parameters::init(argc, argv);

    int failures = 0;

    if (!testSetters())
    {
        std::cout << "Setters of the integration point source failed" << std::endl;
        ++failures;
    }

    if (!testThrowsUnsupportedMethod())
    {
        std::cout << "No exception for an unsupported discretization method" << std::endl;
        ++failures;
    }

    {
        using Grid = Dune::YaspGrid<dim>;
        Grid grid(GlobalPosition(1.0), std::array<int, dim>{{3, 3, 3}});
        failures += testGrid(grid, "cube");

        BoxFVGridGeometry<double, Grid::LeafGridView> gridGeometry(grid.leafGridView());
        if (!testThrowsOutsideElement(gridGeometry))
        {
            std::cout << "No exception for a point source outside of its element" << std::endl;
            ++failures;
        }
    }

#if HAVE_DUNE_ALUGRID
    {
        using Grid = Dune::ALUGrid<dim, dim, Dune::simplex, Dune::conforming>;
        const auto grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid(
            GlobalPosition(0.0), GlobalPosition(1.0), std::array<unsigned int, dim>{{3, 3, 3}}
        );
        failures += testGrid(*grid, "simplex");
    }
#endif

    if (failures > 0)
    {
        std::cout << failures << " checks failed" << std::endl;
        return 1;
    }
    return 0;
}
