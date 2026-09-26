// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test sub-control volume face normals of the control-volume finite element
 *        schemes on a surface grid (dim = 2, dimWorld = 3)
 *
 * On a surface grid the face between two sub-control volumes is a segment, whose
 * two corners do not determine a normal on their own. The normal is fixed by also
 * requiring it to be tangential to the element, which is what is checked here.
 */
#include <config.h>

#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/math.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/facecentered/diamond/fvgridgeometry.hh>
#include <dumux/discretization/pq1bubble/fvgridgeometry.hh>
#include <dumux/discretization/pq2/fvgridgeometry.hh>
#include <dumux/discretization/pq3/fvgridgeometry.hh>

namespace {

using Grid = Dune::FoamGrid<2, 3>;
using GridView = Grid::LeafGridView;
using GlobalPosition = Dune::FieldVector<double, 3>;

constexpr double eps = 1e-12;
int failures = 0;

void check(bool ok, const std::string& message)
{
    if (!ok)
    {
        std::cerr << "FAIL: " << message << std::endl;
        ++failures;
    }
}

std::shared_ptr<Grid> makeGrid(const std::vector<GlobalPosition>& vertices,
                               const std::vector<std::vector<unsigned int>>& elements)
{
    Dune::GridFactory<Grid> factory;
    for (const auto& v : vertices)
        factory.insertVertex(v);
    for (const auto& e : elements)
        factory.insertElement(Dune::GeometryTypes::triangle, e);
    return std::shared_ptr<Grid>(factory.createGrid());
}

//! a single triangle in general position, so every facet is on the domain boundary
std::shared_ptr<Grid> makeTiltedTriangle()
{
    return makeGrid({{0.3, -0.2, 0.5}, {1.7, 0.4, 1.1}, {0.1, 1.3, 2.0}}, {{0, 1, 2}});
}

//! four triangles around an apex, no two of them coplanar
std::shared_ptr<Grid> makeTent()
{
    return makeGrid(
        {{0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.2}, {-1.0, 0.0, 0.1}, {0.0, -1.0, 0.3}},
        {{0, 1, 2}, {0, 2, 3}, {0, 3, 4}, {0, 4, 1}}
    );
}

GlobalPosition elementNormal(const Grid::Codim<0>::Entity::Geometry& geo)
{
    auto n = Dumux::crossProduct(GlobalPosition(geo.corner(1) - geo.corner(0)),
                                 GlobalPosition(geo.corner(2) - geo.corner(0)));
    n /= n.two_norm();
    return n;
}

/*!
 * \brief Check that each normal is a unit vector tangential to its element and
 *        orthogonal to its own face, and points from the inside to the outside
 *        sub-control volume. These four conditions determine the normal uniquely.
 */
template<class GridGeometry>
void checkNormals(const GridView& gridView, const std::string& scheme)
{
    GridGeometry gridGeometry(gridView);
    auto fvGeometry = localView(gridGeometry);

    for (const auto& element : elements(gridView))
    {
        fvGeometry.bind(element);
        const auto en = elementNormal(element.geometry());

        for (const auto& scvf : scvfs(fvGeometry))
        {
            const GlobalPosition n = scvf.unitOuterNormal();
            check(std::abs(n.two_norm() - 1.0) < eps, scheme + ": normal is not a unit vector");
            check(std::abs(n*en) < eps, scheme + ": normal is not tangential to the element");

            const auto faceGeometry = fvGeometry.geometry(scvf);
            const GlobalPosition tangent = faceGeometry.corner(faceGeometry.corners() - 1)
                                         - faceGeometry.corner(0);
            check(std::abs(n*tangent) < eps*tangent.two_norm(),
                  scheme + ": normal is not orthogonal to its face");

            if (!scvf.boundary())
            {
                const GlobalPosition d = fvGeometry.scv(scvf.outsideScvIdx()).dofPosition()
                                       - fvGeometry.scv(scvf.insideScvIdx()).dofPosition();
                check(n*d > 0.0, scheme + ": normal is not oriented from the inside to the outside scv");
            }
        }
    }
}

/*!
 * \brief Check the divergence theorem on each sub-control volume: the integral of the
 *        outer normal over its closed boundary vanishes. This only holds element-wise
 *        when every facet of the element is on the domain boundary, so that the whole
 *        boundary of each sub-control volume is covered by sub-control volume faces.
 */
template<class GridGeometry>
void checkClosure(const GridView& gridView, const std::string& scheme)
{
    GridGeometry gridGeometry(gridView);
    auto fvGeometry = localView(gridGeometry);

    for (const auto& element : elements(gridView))
    {
        fvGeometry.bind(element);
        std::vector<GlobalPosition> outflow(fvGeometry.numScv(), GlobalPosition(0.0));

        for (const auto& scvf : scvfs(fvGeometry))
        {
            GlobalPosition contribution = scvf.unitOuterNormal();
            contribution *= scvf.area();
            outflow[scvf.insideScvIdx()] += contribution;
            if (!scvf.boundary())
                outflow[scvf.outsideScvIdx()] -= contribution;
        }

        for (const auto& o : outflow)
            check(o.two_norm() < eps, scheme + ": sub-control volume boundary does not close ("
                                      + std::to_string(o.two_norm()) + ")");
    }
}

} // end anonymous namespace

int main(int argc, char* argv[])
{
    Dumux::initialize(argc, argv);

    using Box = Dumux::BoxFVGridGeometry<double, GridView, true>;
    using Diamond = Dumux::FaceCenteredDiamondFVGridGeometry<GridView, true>;
    using PQ1Bubble = Dumux::PQ1BubbleFVGridGeometry<double, GridView, true>;
    using HybridPQ1Bubble = Dumux::PQ1BubbleFVGridGeometry<double, GridView, true,
        Dumux::HybridPQ1BubbleCVFEGridGeometryTraits<Dumux::PQ1BubbleDefaultGridGeometryTraits<GridView>>>;
    using PQ2 = Dumux::PQ2FVGridGeometry<double, GridView, true>;
    using PQ3 = Dumux::PQ3FVGridGeometry<double, GridView, true>;

    for (const auto& [grid, name] : {std::pair{makeTiltedTriangle(), std::string{"triangle"}},
                                     std::pair{makeTent(), std::string{"tent"}}})
    {
        const auto gridView = grid->leafGridView();
        checkNormals<Box>(gridView, "box/" + name);
        checkNormals<Diamond>(gridView, "diamond/" + name);
        checkNormals<PQ1Bubble>(gridView, "pq1bubble/" + name);
        checkNormals<HybridPQ1Bubble>(gridView, "hybridpq1bubble/" + name);
        checkNormals<PQ2>(gridView, "pq2/" + name);
        checkNormals<PQ3>(gridView, "pq3/" + name);
    }

    // the sub-control volumes only close on a mesh of one element, where every
    // facet is on the domain boundary
    {
        const auto grid = makeTiltedTriangle();
        const auto gridView = grid->leafGridView();
        checkClosure<Box>(gridView, "box");
        checkClosure<Diamond>(gridView, "diamond");
        checkClosure<HybridPQ1Bubble>(gridView, "hybridpq1bubble");
        checkClosure<PQ2>(gridView, "pq2");
        checkClosure<PQ3>(gridView, "pq3");
    }

    if (failures > 0)
    {
        std::cerr << failures << " check(s) failed" << std::endl;
        return 1;
    }

    std::cout << "All scvf normal checks passed on a dim=2, dimWorld=3 grid" << std::endl;
    return 0;
}
