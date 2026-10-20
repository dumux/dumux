// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test for the limited linear reconstruction of the shallow water state at faces:
 *        exactness for linear fields, bounds on Cartesian grids, independence of the
 *        coordinate system, lake at rest, and the first-order fallback next to dry cells
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/flux/shallowwater/reconstruction.hh>

namespace Dumux::Test {

using Scalar = double;
using GlobalPosition = Dune::FieldVector<Scalar, 2>;
constexpr Scalar gravity = 9.81;

struct State { Scalar h, z, u, v; };

struct VolVars
{
    State s;
    Scalar waterDepth() const { return s.h; }
    Scalar bedSurface() const { return s.z; }
    Scalar velocity(int k) const { return k == 0 ? s.u : s.v; }
};

struct GridVolVars { static constexpr bool cachingEnabled = true; };

struct ElementVolVars
{
    using GridVolumeVariables = GridVolVars;
    std::vector<VolVars> volVars;
    const VolVars& operator[](std::size_t i) const { return volVars[i]; }
};

using Field = std::function<State(const GlobalPosition&)>;

struct FaceResult
{
    std::size_t inside;
    GlobalPosition insideCenter;
    GlobalPosition faceCenter;
    bool boundary;
    ShallowWater::FaceState<Scalar> state;
    std::size_t numNeighbors;
    GlobalPosition unitNormal;
};

template<class GridGeometry>
ElementVolVars cellValues(const GridGeometry& gridGeometry, const Field& field)
{
    ElementVolVars ev;
    ev.volVars.resize(gridGeometry.numScv());
    for (const auto& element : elements(gridGeometry.gridView()))
        ev.volVars[gridGeometry.elementMapper().index(element)].s = field(element.geometry().center());
    return ev;
}

template<class GridGeometry>
std::vector<FaceResult> reconstructAllFaces(const GridGeometry& gridGeometry,
                                            const ElementVolVars& ev,
                                            const ShallowWater::ReconstructionParameters<Scalar>& params)
{
    std::vector<FaceResult> results;
    auto fvGeometry = localView(gridGeometry);
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        fvGeometry.bind(element);
        const auto i = gridGeometry.elementMapper().index(element);
        const ShallowWater::CellReconstruction<Scalar, GlobalPosition> rec(params, fvGeometry, ev, i, gravity);
        std::size_t numNeighbors = 0;
        for (const auto& scvf : scvfs(fvGeometry))
            numNeighbors += !scvf.boundary();
        for (const auto& scvf : scvfs(fvGeometry))
        {
            FaceResult r{i, fvGeometry.scv(i).center(), scvf.center(), scvf.boundary(), {}, numNeighbors, scvf.unitOuterNormal()};
            if (scvf.boundary())
                r.state = rec.faceState(scvf.center(), scvf.unitOuterNormal());
            else
                r.state = rec.faceState(scvf.center(), scvf.unitOuterNormal(),
                                        ev[scvf.outsideScvIdx()], fvGeometry.scv(scvf.outsideScvIdx()).center());
            results.push_back(r);
        }
    }
    return results;
}

ShallowWater::ReconstructionParameters<Scalar> parameters(bool limit, Scalar relativeEpsilon = 1e-3)
{ return {limit, relativeEpsilon, 1e-5, 1e-3}; }

void check(bool condition, const std::string& message)
{
    if (!condition)
        DUNE_THROW(Dune::Exception, message);
}

// triangulation of the unit square (each square cell split along a diagonal) with transformed vertex positions
std::unique_ptr<Dune::UGGrid<2>> makeSimplexGrid(std::array<unsigned int, 2> cells,
                                                 const std::function<GlobalPosition(const GlobalPosition&)>& transform)
{
    using Grid = Dune::UGGrid<2>;
    Dune::GridFactory<Grid> factory;
    const auto [nx, ny] = cells;
    for (unsigned int j = 0; j <= ny; ++j)
        for (unsigned int i = 0; i <= nx; ++i)
            factory.insertVertex(transform(GlobalPosition{Scalar(i)/nx, Scalar(j)/ny}));

    auto vertex = [nx = nx](unsigned int i, unsigned int j) { return j*(nx + 1) + i; };
    for (unsigned int j = 0; j < ny; ++j)
        for (unsigned int i = 0; i < nx; ++i)
        {
            factory.insertElement(Dune::GeometryTypes::triangle, {vertex(i, j), vertex(i+1, j), vertex(i, j+1)});
            factory.insertElement(Dune::GeometryTypes::triangle, {vertex(i+1, j), vertex(i+1, j+1), vertex(i, j+1)});
        }
    return std::unique_ptr<Grid>(factory.createGrid());
}

template<class GridView>
using GridGeometry = CCTpfaFVGridGeometry<GridView, true>;

// the van Albada average equals b for a = b, does not change sign for estimates of opposite signs large compared
// to the smoothing scale, does not limit increments small compared to it, and is continuously differentiable
void checkVanAlbadaAverage()
{
    const Scalar eps2 = 1e-4;
    std::mt19937 gen(7);
    std::uniform_real_distribution<Scalar> dist(-10.0, 10.0);
    for (int i = 0; i < 10000; ++i)
    {
        const Scalar a = dist(gen);
        check(ShallowWater::vanAlbadaAverage(a, a, eps2) == a, "van Albada average is not exact for equal estimates");
        check(ShallowWater::vanAlbadaAverage(a, a, 0.0) == a, "clipped van Albada average is not exact for equal estimates");
        const Scalar b = dist(gen);
        if (a*b < 0.0 && std::abs(a) > 1.0 && std::abs(b) > 1.0)
            check(std::abs(ShallowWater::vanAlbadaAverage(a, b, eps2)) < 1e-3, "increment of opposite estimates is not limited");
        if (a*b < 0.0)
            check(ShallowWater::vanAlbadaAverage(a, b, 0.0) == 0.0, "clipped average of opposite estimates is not zero");
    }

    const Scalar small = 1e-6;
    check(std::abs(ShallowWater::vanAlbadaAverage(small, -0.5*small, eps2) - 0.25*small) < 1e-3*small,
          "increments small compared to the smoothing scale are limited");

    // one-sided difference quotients at a*b = 0 differ by about h times the second derivative (of order
    // 1/eps^2), whereas a kink as in the clipped average gives a difference of the order of b
    for (const Scalar b : {1.0, -0.3, 0.01})
    {
        auto gap = [&](Scalar eps2Test) {
            const Scalar h = 1e-8;
            const Scalar right = (ShallowWater::vanAlbadaAverage(h, b, eps2Test) - ShallowWater::vanAlbadaAverage(0.0, b, eps2Test))/h;
            const Scalar left = (ShallowWater::vanAlbadaAverage(0.0, b, eps2Test) - ShallowWater::vanAlbadaAverage(-h, b, eps2Test))/h;
            return std::abs(right - left);
        };
        check(gap(eps2) < 1e-3, "van Albada average is not differentiable at ab = 0");
        check(gap(0.0) > 0.1*std::abs(b), "clipped van Albada average is expected to have a kink at ab = 0");
    }
    std::cout << "van Albada average: exact for equal estimates, sign-preserving, smooth" << std::endl;
}

// linear fields are reconstructed exactly at all faces of cells with at least two face neighbors, with and
// without limiting, except for the tangential velocity at boundary faces, which is that of the cell; a cell
// with one face neighbor has no gradient information across it
template<class GridView>
void checkLinearExactness(const GridView& gridView, const std::string& name)
{
    GridGeometry<GridView> gridGeometry(gridView);
    const Field field = [](const GlobalPosition& x) -> State {
        return {1.0 + 0.3*x[0] - 0.2*x[1], 0.1 - 0.5*x[0] + 0.7*x[1], 0.4 - 1.1*x[0] + 0.6*x[1], -0.2 + 0.9*x[0] + 0.3*x[1]};
    };
    const auto ev = cellValues(gridGeometry, field);
    for (const bool limit : {false, true})
    {
        Scalar maxError = 0.0, maxErrorOneNeighbor = 0.0;
        for (const auto& r : reconstructAllFaces(gridGeometry, ev, parameters(limit)))
        {
            auto exact = field(r.faceCenter);
            if (r.boundary)
            {
                const auto cell = field(r.insideCenter);
                const Scalar normalIncrement = (exact.u - cell.u)*r.unitNormal[0] + (exact.v - cell.v)*r.unitNormal[1];
                exact.u = cell.u + normalIncrement*r.unitNormal[0];
                exact.v = cell.v + normalIncrement*r.unitNormal[1];
            }
            using std::abs; using std::max;
            auto& m = r.numNeighbors >= 2 ? maxError : maxErrorOneNeighbor;
            m = max({m, abs(r.state.waterDepth - exact.h), abs(r.state.bedSurface - exact.z),
                        abs(r.state.velocityX - exact.u), abs(r.state.velocityY - exact.v)});
        }
        std::cout << name << (limit ? " (limited)" : " (unlimited)") << ": max error of linear fields " << maxError
                  << " (cells with one face neighbor: " << maxErrorOneNeighbor << ")" << std::endl;
        check(maxError < 1e-12, name + ": linear fields are not reconstructed exactly");
    }
}

// on a uniform Cartesian grid, the increment to a face limited without smoothing has the sign of the difference
// to the neighbor across the face and at most (1 + sqrt(2))/4 = 0.604 times its magnitude (the maximum of the
// van Albada average of d^-/2 and d^+/2 over the ratio of the backward and forward differences d^- and d^+)
void checkCartesianBounds()
{
    using Grid = Dune::YaspGrid<2>;
    const auto grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(GlobalPosition(0.0), GlobalPosition(1.0), std::array<unsigned int, 2>{{9, 9}});
    GridGeometry<typename Grid::LeafGridView> gridGeometry(grid->leafGridView());

    std::mt19937 gen(42);
    std::uniform_real_distribution<Scalar> dist(-1.0, 1.0);
    ElementVolVars ev;
    ev.volVars.resize(gridGeometry.numScv());
    for (auto& vv : ev.volVars)
        vv.s = {1.0 + 0.5*dist(gen), dist(gen), dist(gen), dist(gen)};

    std::size_t numChecked = 0;
    auto fvGeometry = localView(gridGeometry);
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        fvGeometry.bind(element);
        const auto i = gridGeometry.elementMapper().index(element);
        const ShallowWater::CellReconstruction<Scalar, GlobalPosition> rec(parameters(true, 0.0), fvGeometry, ev, i, gravity);
        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (scvf.boundary())
                continue;
            const auto j = scvf.outsideScvIdx();
            const auto st = rec.faceState(scvf.center(), scvf.unitOuterNormal(), ev[j], fvGeometry.scv(j).center());
            const auto& si = ev[i].s;
            const auto& sj = ev[j].s;
            auto between = [](Scalar value, Scalar cell, Scalar neighbor) {
                const Scalar bound = cell + 0.61*(neighbor - cell);
                return value >= std::min(cell, bound) - 1e-14 && value <= std::max(cell, bound) + 1e-14;
            };
            check(between(st.waterDepth, si.h, sj.h), "water depth outside of bounds");
            check(between(st.waterDepth + st.bedSurface, si.h + si.z, sj.h + sj.z), "free surface outside of bounds");
            check(between(st.velocityX, si.u, sj.u), "x-velocity outside of bounds");
            check(between(st.velocityY, si.v, sj.v), "y-velocity outside of bounds");
            ++numChecked;
        }
    }
    std::cout << "Cartesian bounds hold at " << numChecked << " interior faces for random data" << std::endl;
}

// the reconstruction on a rotated grid with rotated fields gives the rotated face states
void checkRotation()
{
    const Scalar angle = 0.7;
    const Scalar c = std::cos(angle), s = std::sin(angle);
    auto rotate = [=](const GlobalPosition& x) { return GlobalPosition{c*x[0] - s*x[1], s*x[0] + c*x[1]}; };
    auto rotateBack = [=](const GlobalPosition& x) { return GlobalPosition{c*x[0] + s*x[1], -s*x[0] + c*x[1]}; };
    auto jitter = [](const GlobalPosition& x) {
        return GlobalPosition{x[0] + 0.03*std::sin(17.0*x[1] + 3.0*x[0]), x[1] + 0.03*std::cos(13.0*x[0] - 5.0*x[1])};
    };

    const auto gridA = makeSimplexGrid({{7, 6}}, jitter);
    const auto gridB = makeSimplexGrid({{7, 6}}, [&](const GlobalPosition& x) { return rotate(jitter(x)); });
    GridGeometry<typename Dune::UGGrid<2>::LeafGridView> ggA(gridA->leafGridView()), ggB(gridB->leafGridView());

    const Field fieldA = [](const GlobalPosition& x) -> State {
        return {1.0 + 0.3*std::sin(4.0*x[0])*std::cos(3.0*x[1]), 0.2*std::cos(5.0*x[0] + x[1]),
                std::sin(6.0*x[1]) - x[0]*x[0], std::cos(5.0*x[0])*x[1]};
    };
    const Field fieldB = [&](const GlobalPosition& x) -> State {
        const auto sA = fieldA(rotateBack(x));
        const auto u = rotate(GlobalPosition{sA.u, sA.v});
        return {sA.h, sA.z, u[0], u[1]};
    };

    const auto resultsA = reconstructAllFaces(ggA, cellValues(ggA, fieldA), parameters(true));
    const auto resultsB = reconstructAllFaces(ggB, cellValues(ggB, fieldB), parameters(true));
    check(resultsA.size() == resultsB.size(), "different number of faces");

    Scalar maxError = 0.0;
    for (const auto& a : resultsA)
    {
        const auto insideB = rotate(a.insideCenter);
        const auto faceB = rotate(a.faceCenter);
        auto it = std::find_if(resultsB.begin(), resultsB.end(), [&](const auto& b) {
            return (b.insideCenter - insideB).two_norm() < 1e-10 && (b.faceCenter - faceB).two_norm() < 1e-10;
        });
        check(it != resultsB.end(), "no matching face on the rotated grid");
        const auto u = rotate(GlobalPosition{a.state.velocityX, a.state.velocityY});
        using std::abs; using std::max;
        maxError = max({maxError, abs(it->state.waterDepth - a.state.waterDepth), abs(it->state.bedSurface - a.state.bedSurface),
                                  abs(it->state.velocityX - u[0]), abs(it->state.velocityY - u[1])});
    }
    std::cout << "Rotated grid: max difference of the rotated face states " << maxError << std::endl;
    check(maxError < 1e-12, "reconstruction depends on the coordinate system");
}

// for a constant free surface with a nonlinear bed, the face states keep the free surface and zero velocity
template<class GridView>
void checkLakeAtRest(const GridView& gridView, const std::string& name)
{
    GridGeometry<GridView> gridGeometry(gridView);
    const Field field = [](const GlobalPosition& x) -> State {
        const Scalar z = 0.8*std::exp(-((x[0] - 0.5)*(x[0] - 0.5) + (x[1] - 0.4)*(x[1] - 0.4))/0.05);
        return {2.0 - z, z, 0.0, 0.0};
    };
    Scalar maxError = 0.0;
    for (const auto& r : reconstructAllFaces(gridGeometry, cellValues(gridGeometry, field), parameters(true)))
    {
        using std::abs; using std::max;
        maxError = max({maxError, abs(r.state.waterDepth + r.state.bedSurface - 2.0), abs(r.state.velocityX), abs(r.state.velocityY)});
    }
    std::cout << name << ": max deviation from the lake at rest at faces " << maxError << std::endl;
    check(maxError < 1e-14, name + ": lake at rest is not reconstructed exactly");
}

// cells with a dry cell in their stencil use their cell values at all faces, cells further away do not
void checkDryNeighbor()
{
    using Grid = Dune::YaspGrid<2>;
    const auto grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(GlobalPosition(0.0), GlobalPosition(1.0), std::array<unsigned int, 2>{{7, 7}});
    GridGeometry<typename Grid::LeafGridView> gridGeometry(grid->leafGridView());
    const GlobalPosition dryCenter{0.5, 0.5};
    const Field field = [&](const GlobalPosition& x) -> State {
        if ((x - dryCenter).two_norm() < 1e-8)
            return {0.0, 0.3, 0.0, 0.0};
        return {1.0 + 0.2*x[0], 0.1*x[1], 0.3*x[1], -0.1*x[0]};
    };
    const auto ev = cellValues(gridGeometry, field);
    const Scalar dx = 1.0/7.0;
    for (const auto& r : reconstructAllFaces(gridGeometry, ev, parameters(true)))
    {
        const auto& cell = ev[r.inside].s;
        const bool firstOrder = r.state.waterDepth == cell.h && r.state.bedSurface == cell.z
                                && r.state.velocityX == cell.u && r.state.velocityY == cell.v;
        const bool inStencil = (r.insideCenter - dryCenter).two_norm() < 1.01*dx;
        check(firstOrder == inStencil, "first-order fallback applied to the wrong cells");
    }
    std::cout << "Dry cell: first-order fallback exactly in the cells whose stencil contains the dry cell" << std::endl;
}

} // end namespace Dumux::Test

int main(int argc, char *argv[])
{
    using namespace Dumux;
    using namespace Dumux::Test;
    initialize(argc, argv);

    using Yasp = Dune::YaspGrid<2>;
    const auto yasp = Dune::StructuredGridFactory<Yasp>::createCubeGrid(GlobalPosition(0.0), GlobalPosition{1.0, 0.6}, std::array<unsigned int, 2>{{7, 5}});
    const auto simplex = makeSimplexGrid({{6, 5}}, [](const GlobalPosition& x) { return x; });
    const auto jittered = makeSimplexGrid({{6, 5}}, [](const GlobalPosition& x) {
        return GlobalPosition{x[0] + 0.04*std::sin(11.0*x[1] + 2.0*x[0]), x[1] + 0.04*std::cos(9.0*x[0] - 4.0*x[1])};
    });

    checkVanAlbadaAverage();
    checkLinearExactness(yasp->leafGridView(), "YaspGrid");
    checkLinearExactness(simplex->leafGridView(), "Simplex grid");
    checkLinearExactness(jittered->leafGridView(), "Distorted simplex grid");
    checkCartesianBounds();
    checkRotation();
    checkLakeAtRest(yasp->leafGridView(), "YaspGrid");
    checkLakeAtRest(jittered->leafGridView(), "Distorted simplex grid");
    checkDryNeighbor();

    return 0;
}
