// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test for the trace abstraction of the mortar coupling: the trace concepts, the
 *        trace kind of a subdomain, and the facet trace of bulk subdomains.
 */
#include <config.h>

#include <cmath>
#include <memory>
#include <iostream>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/foamgrid/foamgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dumux/common/initialize.hh>
#include <dumux/common/concepts/entityset_.hh>
#include <dumux/common/concepts/mortartrace_.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/intersectionentityset.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/multidomain/mortar/trace.hh>
#include <dumux/multidomain/mortar/projectors.hh>

namespace {

// a grid-free trace as a network subdomain would provide it: declarations suffice for
// the concept checks
struct WindowLikeTrace
{
    using EntitySet = Dumux::GeometriesEntitySet<Dune::MultiLinearGeometry<double, 1, 2>>;
    std::size_t size() const;
    std::shared_ptr<const EntitySet> entitySet() const;
    const Dumux::BoundingBoxTree<EntitySet>& boundingBoxTree() const;
};

struct TraceWithoutTree
{
    using EntitySet = Dumux::GeometriesEntitySet<Dune::MultiLinearGeometry<double, 1, 2>>;
    std::size_t size() const;
    std::shared_ptr<const EntitySet> entitySet() const;
};

template<int dimWorld>
auto makeMortarGrid(std::size_t numCells, double lowerLeftCoordinate = 0.0, double upperRightCoordinate = 1.0, double height = 0.0)
{
    using Grid = Dune::FoamGrid<dimWorld-1, dimWorld>;
    Dune::GridFactory<Grid> factory;
    if constexpr (dimWorld == 2)
    {
        for (std::size_t i = 0; i <= numCells; ++i)
        {
            const double x = lowerLeftCoordinate + (upperRightCoordinate - lowerLeftCoordinate)*i/numCells;
            factory.insertVertex({x, height});
        }
        for (unsigned int i = 0; i < numCells; ++i)
            factory.insertElement(Dune::GeometryTypes::line, {i, i+1});
    }
    else
    {
        // a two-dimensional foam grid holds triangles only
        for (std::size_t j = 0; j <= numCells; ++j)
            for (std::size_t i = 0; i <= numCells; ++i)
                factory.insertVertex({
                    lowerLeftCoordinate + (upperRightCoordinate - lowerLeftCoordinate)*i/numCells,
                    lowerLeftCoordinate + (upperRightCoordinate - lowerLeftCoordinate)*j/numCells,
                    height
                });
        const auto v = [&] (unsigned int i, unsigned int j) { return static_cast<unsigned int>(j*(numCells+1) + i); };
        for (unsigned int j = 0; j < numCells; ++j)
            for (unsigned int i = 0; i < numCells; ++i)
            {
                factory.insertElement(Dune::GeometryTypes::triangle, {v(i, j), v(i+1, j), v(i+1, j+1)});
                factory.insertElement(Dune::GeometryTypes::triangle, {v(i, j), v(i+1, j+1), v(i, j+1)});
            }
    }
    return std::shared_ptr<Grid>(factory.createGrid());
}

int check(bool condition, const std::string& what)
{
    if (!condition)
        std::cout << "FAILED: " << what << std::endl;
    return condition ? 0 : 1;
}

} // end anonymous namespace

int main(int argc, char** argv)
{
    using namespace Dumux;
    initialize(argc, argv);

    using Grid2 = Dune::YaspGrid<2>;
    using Grid3 = Dune::YaspGrid<3>;
    using MortarGrid2 = Dune::FoamGrid<1, 2>;
    using MortarGrid3 = Dune::FoamGrid<2, 3>;
    using BoxGG2 = BoxFVGridGeometry<double, typename Grid2::LeafGridView>;
    using TpfaGG2 = CCTpfaFVGridGeometry<typename Grid2::LeafGridView>;
    using BoxGG3 = BoxFVGridGeometry<double, typename Grid3::LeafGridView>;
    using MortarGG2 = CCTpfaFVGridGeometry<typename MortarGrid2::LeafGridView>;
    using MortarGG3 = CCTpfaFVGridGeometry<typename MortarGrid3::LeafGridView>;
    using FacetTrace2 = Mortar::FacetTrace<Grid2, MortarGrid2>;
    using FacetTrace3 = Mortar::FacetTrace<Grid3, MortarGrid3>;

    // the entity sets a trace may be built on
    static_assert(Concept::GeometricEntitySet<GridViewGeometricEntitySet<typename MortarGrid2::LeafGridView>>);
    static_assert(Concept::GeometricEntitySet<GeometriesEntitySet<Dune::MultiLinearGeometry<double, 1, 2>>>);
    static_assert(Concept::GeometricEntitySet<GeometriesEntitySet<Dune::MultiLinearGeometry<double, 2, 3>>>);
    static_assert(!Concept::GeometricEntitySet<std::vector<double>>);

    // the trace concepts discriminate: a grid-free trace is a trace but not a grid trace
    static_assert(Concept::MortarTrace<FacetTrace2>);
    static_assert(Concept::GridMortarTrace<FacetTrace2>);
    static_assert(Concept::MortarTrace<FacetTrace3>);
    static_assert(Concept::GridMortarTrace<FacetTrace3>);
    static_assert(Concept::MortarTrace<WindowLikeTrace>);
    static_assert(!Concept::GridMortarTrace<WindowLikeTrace>);
    static_assert(!Concept::MortarTrace<TraceWithoutTree>);
    static_assert(!Concept::MortarTrace<int>);

    // the trace kind of bulk subdomains, and the absence of a kind for a dimension mismatch
    static_assert(Mortar::HasTrace<BoxGG2, MortarGrid2>);
    static_assert(Mortar::HasTrace<TpfaGG2, MortarGrid2>);
    static_assert(Mortar::HasTrace<BoxGG3, MortarGrid3>);
    static_assert(std::is_same_v<Mortar::Trace<BoxGG2, MortarGrid2>, FacetTrace2>);
    static_assert(std::is_same_v<Mortar::Trace<TpfaGG2, MortarGrid2>, FacetTrace2>);
    static_assert(std::is_same_v<Mortar::Trace<BoxGG3, MortarGrid3>, FacetTrace3>);
    static_assert(!Mortar::HasTrace<BoxGG2, MortarGrid3>);
    static_assert(!Mortar::HasTrace<BoxFVGridGeometry<double, typename MortarGrid2::LeafGridView>, MortarGrid2>);

    int exitCode = 0;

    // a bulk subdomain in 2d against a mortar of three cells along its lower boundary
    Grid2 grid2{{1.0, 1.0}, {4, 4}};
    auto gg2 = std::make_shared<BoxGG2>(grid2.leafGridView());
    auto mortarGrid2 = makeMortarGrid<2>(3);
    auto mortarGG2 = std::make_shared<MortarGG2>(mortarGrid2->leafGridView());

    const FacetTrace2 trace2{*gg2, *mortarGG2};
    exitCode += check(trace2.size() == 4, "2d facet trace has the four lower boundary facets");
    exitCode += check(trace2.entitySet()->size() == trace2.size(), "entity set size equals trace size");
    exitCode += check(&trace2.boundingBoxTree().entitySet() == trace2.entitySet().get(), "bounding box tree is built over the trace's entity set");
    exitCode += check(&trace2.boundingBoxTree() == &trace2.boundingBoxTree(), "bounding box tree is built once");
    exitCode += check(trace2.gridView().size(0) == 4, "grid view exposes the trace cells");
    {
        bool onLowerBoundary = true;
        bool consistentIndices = true;
        bool oneHostIntersectionEach = true;
        std::vector<bool> seen(trace2.size(), false);
        for (const auto& cell : *trace2.entitySet())
        {
            onLowerBoundary &= std::abs(cell.geometry().center()[1]) < 1e-12;
            const auto index = trace2.entitySet()->index(cell);
            consistentIndices &= index < trace2.size() && !seen[index];
            seen[index] = true;
            consistentIndices &= trace2.entitySet()->index(trace2.entitySet()->entity(index)) == index;
            oneHostIntersectionEach &= trace2.hostGridIntersections(cell).size() == 1;
        }
        exitCode += check(onLowerBoundary, "every trace cell lies on the lower boundary");
        exitCode += check(consistentIndices, "trace cells are indexed contiguously by the entity set");
        exitCode += check(oneHostIntersectionEach, "a boundary trace cell has one host intersection");
    }

    // a mortar along an interior grid line is met by no boundary facet
    auto interiorMortarGrid = makeMortarGrid<2>(3, 0.0, 1.0, 0.5);
    auto interiorMortarGG = std::make_shared<MortarGG2>(interiorMortarGrid->leafGridView());
    const FacetTrace2 interiorTrace{*gg2, *interiorMortarGG};
    exitCode += check(interiorTrace.size() == 0, "a mortar on an interior grid line has an empty facet trace");

    // a mortar shorter than the boundary selects only the facets it meets
    auto shortMortarGrid = makeMortarGrid<2>(2, 0.0, 0.5);
    auto shortMortarGG = std::make_shared<MortarGG2>(shortMortarGrid->leafGridView());
    const FacetTrace2 shortTrace{*gg2, *shortMortarGG};
    exitCode += check(shortTrace.size() == 2, "a mortar over half the boundary selects two facets");

    // the selector constructor is not bound to a mortar
    const FacetTrace2 lShaped{grid2, [] (const auto&, const auto& is) {
        const auto c = is.geometry().center();
        return c[1] < 1e-10 || c[0] < 1e-10;
    }};
    exitCode += check(lShaped.size() == 8, "selector constructor picks the eight facets of an L-shaped trace");

    // the glue between mortar and trace, built through the trace interface alone, covers the mortar
    {
        using MortarEntitySet = std::remove_cvref_t<decltype(mortarGG2->boundingBoxTree())>::EntitySet;
        IntersectionEntitySet<MortarEntitySet, FacetTrace2::EntitySet> glue;
        glue.build(mortarGG2->boundingBoxTree(), trace2.boundingBoxTree());
        double measure = 0.0;
        for (const auto& is : intersections(glue))
            measure += is.geometry().volume();
        exitCode += check(glue.size() == 6, "three mortar cells against four trace cells give six intersections");
        exitCode += check(std::abs(measure - 1.0) < 1e-12, "the intersections tile the mortar");
    }

    // trace spaces from the interface: piecewise constant on any trace, continuous on a grid trace
    {
        const auto p0 = Mortar::Detail::makeTraceSpace<0>(trace2);
        const auto p1 = Mortar::Detail::makeTraceSpace<1>(trace2);
        exitCode += check(p0.basis().size() == 4, "piecewise constant trace space has one dof per trace cell");
        exitCode += check(p1.basis().size() == 5, "continuous trace space has one dof per trace vertex");
    }

    // the two-dimensional foam grid holds triangles only, so a hexahedral host grid has no
    // facet trace on it, and the construction says so rather than dropping corners
    auto mortarGrid3 = makeMortarGrid<3>(1);
    auto mortarGG3 = std::make_shared<MortarGG3>(mortarGrid3->leafGridView());
    {
        Grid3 grid3{{1.0, 1.0, 1.0}, {2, 2, 2}};
        auto gg3 = std::make_shared<BoxGG3>(grid3.leafGridView());
        bool threw = false;
        try { const FacetTrace3 trace3{*gg3, *mortarGG3}; }
        catch (const Dune::GridError&) { threw = true; }
        exitCode += check(threw, "quadrilateral facets on a triangle facet grid are refused");
    }

#if HAVE_DUNE_ALUGRID
    // a simplex bulk subdomain in 3d against a mortar plane of two triangles
    {
        using SimplexGrid3 = Dune::ALUGrid<3, 3, Dune::simplex, Dune::nonconforming>;
        using SimplexGG3 = BoxFVGridGeometry<double, typename SimplexGrid3::LeafGridView>;
        using SimplexTrace3 = Mortar::FacetTrace<SimplexGrid3, MortarGrid3>;
        static_assert(std::is_same_v<Mortar::Trace<SimplexGG3, MortarGrid3>, SimplexTrace3>);

        auto grid3 = Dune::StructuredGridFactory<SimplexGrid3>::createSimplexGrid({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, {2, 2, 2});
        auto gg3 = std::make_shared<SimplexGG3>(grid3->leafGridView());
        const SimplexTrace3 trace3{*gg3, *mortarGG3};
        exitCode += check(trace3.size() == 8, "3d facet trace has the eight lower boundary triangles");
        bool onLowerBoundary = true;
        double area = 0.0;
        for (const auto& cell : *trace3.entitySet())
        {
            onLowerBoundary &= std::abs(cell.geometry().center()[2]) < 1e-12;
            area += cell.geometry().volume();
        }
        exitCode += check(onLowerBoundary, "every 3d trace cell lies on the mortar plane");
        exitCode += check(std::abs(area - 1.0) < 1e-12, "the 3d trace cells cover the mortar plane");
        exitCode += check(Mortar::Detail::makeTraceSpace<0>(trace3).basis().size() == 8, "3d piecewise constant trace space has one dof per cell");
        exitCode += check(Mortar::Detail::makeTraceSpace<1>(trace3).basis().size() == 9, "3d continuous trace space has nine vertex dofs");

        using MortarEntitySet = std::remove_cvref_t<decltype(mortarGG3->boundingBoxTree())>::EntitySet;
        IntersectionEntitySet<MortarEntitySet, SimplexTrace3::EntitySet> glue;
        glue.build(mortarGG3->boundingBoxTree(), trace3.boundingBoxTree());
        double measure = 0.0;
        for (const auto& is : intersections(glue))
            measure += is.geometry().volume();
        exitCode += check(std::abs(measure - 1.0) < 1e-12, "the 3d intersections tile the mortar plane");
    }
#endif

    if (exitCode == 0)
        std::cout << "All trace checks passed" << std::endl;
    return exitCode;
}
