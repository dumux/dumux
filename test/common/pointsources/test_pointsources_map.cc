// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test the point source maps for control-volume, hybrid and finite element discretizations.
 *
 * Points in the interior of elements and on their faces, edges and vertices are registered.
 * Each point has to be counted once in total on the control volumes (if the discretization has
 * any) and once in total on the elements. The positions in local element coordinates have to map
 * back to the point. Contexts that know their element have to give the same maps as the search
 * of the grid.
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/pointsources.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/pq1.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/pq2.hh>
#include <dumux/geometry/intersectingentities.hh>
#include <dumux/geometry/intersectspointgeometry.hh>

namespace Dumux {

namespace Properties::TTag {
template<class Discretization, class GridType>
struct PointSourceMapTest { using InheritsFrom = std::tuple<Discretization>; };
} // end namespace Properties::TTag

namespace Properties {

template<class TypeTag, class Discretization, class GridType>
struct Grid<TypeTag, TTag::PointSourceMapTest<Discretization, GridType>>
{ using type = GridType; };

template<class TypeTag, class Discretization, class GridType>
struct Scalar<TypeTag, TTag::PointSourceMapTest<Discretization, GridType>>
{ using type = double; };

template<class TypeTag, class Discretization, class GridType>
struct EnableGridGeometryCache<TypeTag, TTag::PointSourceMapTest<Discretization, GridType>>
{ static constexpr bool value = true; };

} // end namespace Properties

namespace Test {

struct Checker
{
    int failures = 0;

    void check(bool ok, const std::string& what)
    {
        if (!ok)
        {
            ++failures;
            std::cerr << "FAILED: " << what << std::endl;
        }
    }
};

template<class Grid>
std::shared_ptr<Grid> makeGrid()
{
    constexpr int dim = Grid::dimension;
    Dune::FieldVector<double, dim> lower(0.0), upper(1.0);
    upper[0] = 1.5;
    std::array<unsigned int, dim> cells;
    cells.fill(dim == 2 ? 3 : 2);
    cells[0] += 1;
#if HAVE_DUNE_UGGRID
    if constexpr (std::is_same_v<Grid, Dune::UGGrid<dim>>)
        return Dune::StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, cells);
    else
#endif
        return Dune::StructuredGridFactory<Grid>::createCubeGrid(lower, upper, cells);
}

//! a random point in the relative interior of the sub-entity with the given index, in the local coordinates of the element
template<int codim, class Element, class Rng>
auto randomPointOnSubEntity(const Element& element, int subEntity, Rng& rng)
{
    const auto referenceElement = Dune::referenceElement(element.geometry());
    const auto subGeometry = referenceElement.template geometry<codim>(subEntity);
    constexpr int subDim = Element::mydimension - codim;
    Dune::FieldVector<double, subDim> local(0.0);
    std::uniform_real_distribution<double> uniform(0.1, 0.9);
    if constexpr (subDim > 0)
    {
        if (subGeometry.type().isSimplex())
        {
            std::array<double, subDim + 1> weights;
            for (auto& w : weights)
                w = -std::log(uniform(rng));
            const auto sum = std::accumulate(weights.begin(), weights.end(), 0.0);
            for (int i = 0; i < subDim; ++i)
                local[i] = weights[i]/sum;
        }
        else
            for (int i = 0; i < subDim; ++i)
                local[i] = uniform(rng);
    }
    return subGeometry.global(local);
}

template<class Element>
bool isInsideReferenceElement(const Element& element, const typename Element::Geometry::LocalCoordinate& local)
{
    constexpr double eps = 1e-10;
    double sum = 0.0;
    for (const auto& x : local)
    {
        if (x < -eps)
            return false;
        sum += x;
    }
    if (element.type().isSimplex())
        return sum < 1.0 + eps;
    return std::all_of(local.begin(), local.end(), [&](double x) { return x < 1.0 + eps; });
}

//! compare two maps with the same keys entry by entry
template<class MapA, class MapB, class EntryEqual>
bool mapsEqual(const MapA& a, const MapB& b, EntryEqual&& entryEqual)
{
    if (a.size() != b.size())
        return false;
    auto itA = a.begin();
    auto itB = b.begin();
    for (; itA != a.end(); ++itA, ++itB)
    {
        if (itA->first != itB->first || itA->second.size() != itB->second.size())
            return false;
        for (std::size_t i = 0; i < itA->second.size(); ++i)
            if (!entryEqual(itA->second[i], itB->second[i]))
                return false;
    }
    return true;
}

template<class TypeTag>
void testMaps(Checker& checker, const std::string& name, std::mt19937& rng)
{
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridDiscretization::LocalView;
    using GlobalPosition = typename GridDiscretization::GlobalCoordinate;
    constexpr int dim = Grid::dimension;
    constexpr bool hasScvs = Experimental::Concepts::FVElementDiscretization<ElementDiscretization>;
    constexpr bool hasElementContexts = Dumux::Detail::LocalDofs::hasNonCVLocalDofsInterface<ElementDiscretization>();

    auto grid = makeGrid<Grid>();
    const GridDiscretization gridDiscretization(grid->leafGridView());
    const auto& gridView = gridDiscretization.gridView();

    using Base = Experimental::PointSourceContext<GlobalPosition>;
    using Context = Experimental::IdPointSourceContext<Base, std::size_t>;
    using LocatedContext = Experimental::ElementIndexPointSourceContext<Context>;

    // points in the interior (codim 0), on faces (codim 1), ..., on vertices (codim dim)
    std::vector<Context> contexts;
    std::vector<LocatedContext> interiorContexts;
    std::uniform_int_distribution<std::size_t> randomElement(0, gridView.size(0) - 1);
    for (int codim = 0; codim <= dim; ++codim)
        for (int n = 0; n < 12; ++n)
        {
            const auto element = gridDiscretization.element(randomElement(rng));
            const auto referenceElement = Dune::referenceElement(element.geometry());
            const int subEntity = std::uniform_int_distribution<int>(0, referenceElement.size(codim) - 1)(rng);
            typename Grid::template Codim<0>::Entity::Geometry::LocalCoordinate local;
            Dune::Hybrid::forEach(std::make_index_sequence<dim + 1>{}, [&](auto c)
            {
                if (int(c) == codim)
                    local = randomPointOnSubEntity<int(c)>(element, subEntity, rng);
            });

            Context context(Base(element.geometry().global(local)), contexts.size());
            contexts.push_back(context);
            if (codim == 0)
                interiorContexts.emplace_back(context, gridDiscretization.elementMapper().index(element));
        }

    const auto scvMap = Experimental::makePointSourceMap(gridDiscretization, contexts);
    const auto elementMap = Experimental::makeElementPointSourceMap(gridDiscretization, contexts);

    // every point is counted once in total on the control volumes and once on the elements
    std::vector<double> scvWeights(contexts.size(), 0.0), elementWeights(contexts.size(), 0.0);
    for (const auto& [key, scvContexts] : scvMap)
        for (const auto& context : scvContexts)
        {
            scvWeights[context.id()] += 1.0/context.embeddings();
            checker.check(context.scvIndex() == key.second, name + ": scv index of the context is its key");
        }

    for (const auto& [eIdx, elementContexts] : elementMap)
    {
        const auto element = gridDiscretization.element(eIdx);
        for (const auto& context : elementContexts)
        {
            elementWeights[context.id()] += 1.0/context.embeddings();

            const auto global = element.geometry().global(context.localPosition());
            auto difference = global - context.position();
            checker.check(difference.two_norm() < 1e-13*(1.0 + context.position().two_norm()),
                          name + ": local position maps back to the point");
            checker.check(isInsideReferenceElement(element, context.localPosition()),
                          name + ": local position inside the reference element");
        }
    }

    for (std::size_t i = 0; i < contexts.size(); ++i)
    {
        checker.check(std::abs(elementWeights[i] - (hasElementContexts ? 1.0 : 0.0)) < 1e-14,
                      name + ": element weights of point " + std::to_string(i) + " sum to " + std::to_string(elementWeights[i]));
        checker.check(std::abs(scvWeights[i] - (hasScvs ? 1.0 : 0.0)) < 1e-14,
                      name + ": scv weights of point " + std::to_string(i) + " sum to " + std::to_string(scvWeights[i]));
    }

    if constexpr (!hasScvs)
        checker.check(scvMap.empty(), name + ": no scv contexts without scvs");

    // contexts that know their element
    const auto sameScvEntry = [](const auto& a, const auto& b)
    { return a.id() == b.id() && a.embeddings() == b.embeddings() && a.scvIndex() == b.scvIndex(); };
    const auto sameElementEntry = [](const auto& a, const auto& b)
    {
        auto difference = a.localPosition() - b.localPosition();
        return a.id() == b.id() && a.embeddings() == b.embeddings() && difference.two_norm() < 1e-14;
    };

    // points that lie in the interior of a known element
    std::vector<Context> interiorPlain(interiorContexts.begin(), interiorContexts.end());
    checker.check(mapsEqual(Experimental::makePointSourceMap(gridDiscretization, interiorPlain),
                            Experimental::makePointSourceMap(gridDiscretization, interiorContexts), sameScvEntry),
                  name + ": scv map of contexts with known element, interior points");
    checker.check(mapsEqual(Experimental::makeElementPointSourceMap(gridDiscretization, interiorPlain),
                            Experimental::makeElementPointSourceMap(gridDiscretization, interiorContexts), sameElementEntry),
                  name + ": element map of contexts with known element, interior points");

    // points shared by several elements, with the elements found by the search
    std::vector<LocatedContext> locatedContexts;
    for (const auto& context : contexts)
    {
        const auto elements = intersectingEntities(context.position(), gridDiscretization.boundingBoxTree());
        for (const auto eIdx : elements)
        {
            auto copy = context;
            copy.setEmbeddings(elements.size());
            locatedContexts.emplace_back(copy, eIdx);
        }
    }
    checker.check(mapsEqual(scvMap, Experimental::makePointSourceMap(gridDiscretization, locatedContexts), sameScvEntry),
                  name + ": scv map of contexts with known element, all points");
    checker.check(mapsEqual(elementMap, Experimental::makeElementPointSourceMap(gridDiscretization, locatedContexts), sameElementEntry),
                  name + ": element map of contexts with known element, all points");

    // a point outside of the grid is ignored
    GlobalPosition outside(10.0);
    const std::vector<Context> outsideContexts{Context(Base(outside), 0)};
    checker.check(Experimental::makePointSourceMap(gridDiscretization, outsideContexts).empty()
                  && Experimental::makeElementPointSourceMap(gridDiscretization, outsideContexts).empty(),
                  name + ": no contexts for a point outside of the grid");

    std::cout << name << ": " << contexts.size() << " points, " << scvMap.size() << " scv entries, "
              << elementMap.size() << " elements, " << locatedContexts.size() << " located contexts" << std::endl;
}

template<class Discretization>
void runGrids(const std::string& name, Checker& checker, std::mt19937& rng)
{
    using TTag = Properties::TTag::PointSourceMapTest<Discretization, Dune::YaspGrid<2>>;
    testMaps<TTag>(checker, name + " cube 2D", rng);
    testMaps<Properties::TTag::PointSourceMapTest<Discretization, Dune::YaspGrid<3>>>(checker, name + " cube 3D", rng);
#if HAVE_DUNE_UGGRID
    testMaps<Properties::TTag::PointSourceMapTest<Discretization, Dune::UGGrid<2>>>(checker, name + " simplex 2D", rng);
    testMaps<Properties::TTag::PointSourceMapTest<Discretization, Dune::UGGrid<3>>>(checker, name + " simplex 3D", rng);
#endif
}

} // end namespace Test
} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;
    using namespace Dumux::Properties;

    initialize(argc, argv);

    Test::Checker checker;
    std::mt19937 rng(20261113);

    Test::runGrids<TTag::BoxModel>("box", checker, rng);
    Test::runGrids<TTag::PQ1BubbleHybridModel>("pq1bubble hybrid", checker, rng);
    Test::runGrids<TTag::PQ2HybridModel>("pq2 hybrid", checker, rng);
    Test::runGrids<TTag::PQ1FEModel>("pq1 fe", checker, rng);
    Test::runGrids<TTag::PQ2FEModel>("pq2 fe", checker, rng);

    if (checker.failures > 0)
    {
        std::cerr << checker.failures << " checks failed" << std::endl;
        return 1;
    }

    std::cout << "All checks passed" << std::endl;
    return 0;
}
