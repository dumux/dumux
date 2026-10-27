// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test the hybrid control-volume finite element geometries of order two and three
 *        on line grids, in one dimension and embedded in three dimensions
 *
 * Every vertex owns the half elements next to it, so each element has one interior
 * sub-control volume face, a point at its midpoint, whose normal is the unit tangent of
 * the element. The remaining degrees of freedom of an element are its interior Lagrange
 * nodes, which belong to that element only.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/pq2/fvgridgeometry.hh>
#include <dumux/discretization/pq3/fvgridgeometry.hh>

namespace {

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

using FoamGrid = Dune::FoamGrid<1, 3>;
using FoamPosition = Dune::FieldVector<double, 3>;

std::shared_ptr<FoamGrid> makeFoamGrid(const std::vector<FoamPosition>& vertices,
                                       const std::vector<std::vector<unsigned int>>& elements)
{
    Dune::GridFactory<FoamGrid> factory;
    for (const auto& v : vertices)
        factory.insertVertex(v);
    for (const auto& e : elements)
        factory.insertElement(Dune::GeometryTypes::line, e);
    return std::shared_ptr<FoamGrid>(factory.createGrid());
}

std::shared_ptr<FoamGrid> makePolyline()
{
    return makeFoamGrid({{0.0, 0.0, 0.0}, {0.4, 0.1, 0.2}, {0.7, 0.5, 0.1}, {0.8, 0.9, -0.3}},
                        {{0, 1}, {1, 2}, {3, 2}});
}

std::shared_ptr<FoamGrid> makeJunction()
{
    return makeFoamGrid(
        {{0.0, 0.0, 0.0}, {0.5, 0.0, 0.0}, {1.0, 0.0, 0.0},
         {-0.25, 0.4, 0.05}, {-0.5, 0.8, 0.1}, {-0.15, -0.3, 0.35}, {-0.3, -0.6, 0.7}},
        {{0, 1}, {1, 2}, {0, 3}, {3, 4}, {5, 0}, {6, 5}}
    );
}

std::shared_ptr<Dune::OneDGrid> makeOneDGrid()
{ return std::make_shared<Dune::OneDGrid>(std::vector<double>{0.0, 0.3, 0.45, 1.0}); }

/*!
 * \brief Check the control volumes, faces, degrees of freedom and their positions
 * \param order the polynomial order of the scheme (1 for box)
 */
template<class GridGeometry, class GridView>
void checkLineGeometry(const GridView& gridView, int order, const std::string& name)
{
    using GlobalPosition = typename GridView::template Codim<0>::Geometry::GlobalCoordinate;

    GridGeometry gridGeometry(gridView);
    auto fvGeometry = localView(gridGeometry);

    // number of elements at each vertex, to identify the free ends of the line grid
    std::vector<int> valence(gridView.size(GridView::dimension), 0);
    for (const auto& element : elements(gridView))
        for (int i = 0; i < 2; ++i)
            ++valence[gridView.indexSet().subIndex(element, i, GridView::dimension)];

    std::map<std::size_t, GlobalPosition> dofPositions;
    for (const auto& element : elements(gridView))
    {
        fvGeometry.bind(element);
        const auto geometry = element.geometry();
        const GlobalPosition a = geometry.corner(0);
        const GlobalPosition b = geometry.corner(1);
        GlobalPosition tangent = b - a;
        const double length = tangent.two_norm();
        tangent /= length;

        check(fvGeometry.numScv() == 2, name + ": expected two sub-control volumes per element");
        double volume = 0.0;
        for (const auto& scv : scvs(fvGeometry))
            volume += scv.volume();
        check(std::abs(volume - length) < eps*length, name + ": sub-control volumes do not fill the element");

        std::size_t numInterior = 0;
        for (const auto& scvf : scvfs(fvGeometry))
        {
            const GlobalPosition n = scvf.unitOuterNormal();
            check(std::abs(n.two_norm() - 1.0) < eps, name + ": normal is not a unit vector");
            check(std::abs(std::abs(n*tangent) - 1.0) < eps, name + ": normal is not tangential to the element");
            check(std::abs(scvf.area() - 1.0) < eps, name + ": a point face has area one");

            const GlobalPosition inside = fvGeometry.scv(scvf.insideScvIdx()).dofPosition();
            if (!scvf.boundary())
            {
                ++numInterior;
                const GlobalPosition outside = fvGeometry.scv(scvf.outsideScvIdx()).dofPosition();
                check(n*(outside - inside) > 0.0, name + ": normal does not point from the inside to the outside scv");
                const GlobalPosition midpoint = geometry.center();
                check((scvf.ipGlobal() - midpoint).two_norm() < eps*length, name + ": interior face is not at the element midpoint");
            }
            else
            {
                const auto vIdx = gridView.indexSet().subIndex(element, scvf.insideScvIdx(), GridView::dimension);
                check(valence[vIdx] == 1, name + ": boundary face at a vertex that is not a free end");
                check((scvf.ipGlobal() - inside).two_norm() < eps*length, name + ": boundary face is not at its vertex");
                const GlobalPosition other = geometry.center();
                check(n*(inside - other) > 0.0, name + ": boundary normal does not point out of the element");
            }
        }
        check(numInterior == 1, name + ": expected one interior face per element");

        // the interior Lagrange nodes of the element, ordered along it
        std::vector<double> nodes;
        std::size_t numNonCV = 0;
        if constexpr (requires { nonCVLocalDofs(fvGeometry); })
        {
            for (const auto& localDof : nonCVLocalDofs(fvGeometry))
            {
                ++numNonCV;
                const GlobalPosition x = ipData(fvGeometry, localDof).global();
                const double s = (x - a)*tangent;
                check(((x - a) - s*tangent).two_norm() < eps*length, name + ": interior node is not on the element");
                nodes.push_back(s/length);
            }
        }
        check(numNonCV == std::size_t(order - 1), name + ": wrong number of interior degrees of freedom");
        std::sort(nodes.begin(), nodes.end());
        for (std::size_t j = 0; j < nodes.size(); ++j)
            check(std::abs(nodes[j] - double(j + 1)/order) < eps, name + ": interior node not equidistant");

        // each local basis function is one at its own node and zero at the others
        const auto& localBasis = fvGeometry.feLocalBasis();
        std::vector<Dune::FieldVector<double, 1>> values;
        for (const auto& localDof : localDofs(fvGeometry))
        {
            const GlobalPosition x = ipData(fvGeometry, localDof).global();
            localBasis.evaluateFunction(geometry.local(x), values);
            for (const auto& other : localDofs(fvGeometry))
            {
                const double expected = (other.index() == localDof.index()) ? 1.0 : 0.0;
                check(std::abs(values[other.index()] - expected) < eps,
                      name + ": basis function is not nodal at its degree of freedom");
            }

            const auto dofIdx = localDof.dofIndex();
            const auto [it, inserted] = dofPositions.emplace(dofIdx, x);
            if (!inserted)
                check((it->second - x).two_norm() < eps, name + ": one global dof at two positions");
        }
    }

    const std::size_t expected = gridView.size(GridView::dimension) + (order - 1)*gridView.size(0);
    check(gridGeometry.numDofs() == expected, name + ": wrong number of global degrees of freedom");
    check(dofPositions.size() == expected, name + ": not every global degree of freedom is reached");
}

template<class Grid>
void checkAllSchemes(const Grid& grid, const std::string& name)
{
    using GridView = typename Grid::LeafGridView;
    const auto gridView = grid.leafGridView();
    checkLineGeometry<Dumux::BoxFVGridGeometry<double, GridView, true>>(gridView, 1, "box " + name);
    checkLineGeometry<Dumux::PQ2FVGridGeometry<double, GridView, true>>(gridView, 2, "pq2 " + name);
    checkLineGeometry<Dumux::PQ3FVGridGeometry<double, GridView, true>>(gridView, 3, "pq3 " + name);
    checkLineGeometry<Dumux::PQ2FVGridGeometry<double, GridView, false>>(gridView, 2, "pq2 (no cache) " + name);
    checkLineGeometry<Dumux::PQ3FVGridGeometry<double, GridView, false>>(gridView, 3, "pq3 (no cache) " + name);
}

} // end namespace

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);

    checkAllSchemes(*makeOneDGrid(), "interval");
    checkAllSchemes(*makePolyline(), "polyline in 3d");
    checkAllSchemes(*makeJunction(), "junction in 3d");

    if (failures > 0)
    {
        std::cerr << failures << " checks failed" << std::endl;
        return 1;
    }
    std::cout << "all checks passed" << std::endl;
    return 0;
}
