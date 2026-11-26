// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief The grids, the distribution of the network, the pressure fields and the keys of the degrees of freedom
 *        of the tests of the embedded coupling on distributed grids
 */
#ifndef DUMUX_TEST_EMBEDDED_1D3D_PARALLEL_SETUP_HH
#define DUMUX_TEST_EMBEDDED_1D3D_PARALLEL_SETUP_HH

#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/foamgrid/parallel/partition.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/multidomain/embedded/circlepoints.hh>
#include <dumux/multidomain/embedded/lowdimpartition.hh>

namespace Dumux {

//! a position rounded so that it identifies a degree of freedom in the sequential and the distributed grid
using DofKey = std::array<long long, 3>;

template<class Position>
DofKey dofKey(const Position& x)
{ return {std::llround(x[0]*1e9), std::llround(x[1]*1e9), std::llround(x[2]*1e9)}; }

template<class Position>
double bulkPressure(const Position& x, int state)
{ return 1.0 + 0.5*std::sin(2.0*x[0] + state)*std::cos(3.0*x[1]) + 0.25*x[2]*(1.0 + state); }

template<class Position>
double lowDimPressure(const Position& x, int state)
{ return 2.0 + 0.3*std::cos(5.0*x[0] + x[1] - state) - 0.2*x[2]; }

/*!
 * \brief A network of a helix, a branch leaving the unit cube, a straight vessel at generic positions,
 *        a straight vessel on faces and one on edges of the tissue cells
 *
 * Optionally, vessels parallel to the y-axis have quadrature points or points of the circle around them
 * just outside of faces of the tissue cells, within the tolerance of the search for the elements of a point.
 */
template<class Grid>
std::unique_ptr<Grid> makeNetwork()
{
    using Position = Dune::FieldVector<double, 3>;
    Dune::GridFactory<Grid> factory;
    unsigned int numVertices = 0;
    const auto addPolyline = [&](const std::vector<Position>& points, int startVertex)
    {
        int previous = startVertex;
        for (std::size_t i = 0; i < points.size(); ++i)
        {
            factory.insertVertex(points[i]);
            const int current = numVertices++;
            if (previous >= 0)
                factory.insertElement(Dune::GeometryTypes::line, {static_cast<unsigned int>(previous), static_cast<unsigned int>(current)});
            previous = current;
        }
    };

    const int numHelixSegments = 40;
    std::vector<Position> helix;
    for (int i = 0; i <= numHelixSegments; ++i)
    {
        const double t = double(i)/numHelixSegments;
        const double phi = 0.37 + 2.0*M_PI*1.3*t;
        helix.push_back({0.5 + 0.23*std::cos(phi), 0.5 + 0.23*std::sin(phi), 0.06 + 0.88*t});
    }
    addPolyline(helix, -1);

    const auto lerp = [](const Position& a, const Position& b, double t){ auto x = a; x.axpy(t, b - a); return x; };
    std::vector<Position> branch;
    const Position branchEnd{0.93, 0.11, 0.37}, outside{1.07, 0.04, 0.29};
    for (int i = 1; i <= 8; ++i)
        branch.push_back(lerp(helix[20], branchEnd, i/8.0));
    for (int i = 1; i <= 2; ++i)
        branch.push_back(lerp(branchEnd, outside, i/2.0));
    addPolyline(branch, 20);

    std::vector<Position> straight;
    const Position a{0.13, 0.81, 0.07}, b{0.71, 0.27, 0.88};
    for (int i = 0; i <= 15; ++i)
        straight.push_back(lerp(a, b, i/15.0));
    addPolyline(straight, -1);

    // a vessel in the plane x = 0.5, which contains faces of the tissue cells, so its pieces lie in two cells
    std::vector<Position> onFaces;
    const Position c{0.5, 0.17, 0.23}, d{0.5, 0.83, 0.71};
    for (int i = 0; i <= 11; ++i)
        onFaces.push_back(lerp(c, d, i/11.0));
    addPolyline(onFaces, -1);

    // a vessel on the edge x = y = 0.5 of the tissue cells, so its pieces lie in four cells
    std::vector<Position> onEdges;
    const Position e{0.5, 0.5, 0.13}, f{0.5, 0.5, 0.87};
    for (int i = 0; i <= 9; ++i)
        onEdges.push_back(lerp(e, f, i/9.0));
    addPolyline(onEdges, -1);

    if (getParam<bool>("Network.NearFaceVessels", false))
    {
        const auto cells = getParamFromGroup<std::array<int, 3>>("Tissue", "Grid.Cells");
        const auto numCirclePoints = getParam<int>("MixedDimension.NumCircleSegments");
        const auto circle = EmbeddedCoupling::circlePoints(Position(0.0), Position{0.0, 1.0, 0.0}, 1.0, numCirclePoints);
        using SpatialParams = GetPropType<Properties::TTag::LOWDIMTYPETAG, Properties::SpatialParams>;
        for (const int normal : {0, 2})
        {
            const int other = 2 - normal;
            const auto circleOffset = std::max_element(circle.begin(), circle.end(), [&](const auto& x, const auto& y)
            { return std::abs(x[normal]) < std::abs(y[normal]); })->operator[](normal);

            int count = 0;
            for (const double sign : {-1.0, 1.0})
            {
                for (const double fraction : {0.3, 0.6, 0.9})
                {
                    const auto facePosition = 0.5 + sign*fraction*1e-7/cells[normal];
                    for (const bool circleNearFace : {false, true})
                    {
                        Position p;
                        p[other] = (circleNearFace ? 0.633 : 0.133) + 0.0317*count;
                        p[normal] = facePosition;
                        if (circleNearFace)
                            for (int i = 0; i < 100; ++i)
                                p[normal] = facePosition - SpatialParams::radiusAtPos(p)*circleOffset;

                        std::vector<Position> vessel;
                        for (int i = 0; i <= 6; ++i)
                        {
                            p[1] = 0.17 + 0.11*i;
                            vessel.push_back(p);
                        }
                        addPolyline(vessel, -1);
                    }
                    ++count;
                }
            }
        }
    }

    return factory.createGrid();
}

template<class Grid>
inline constexpr bool isTensorProductYaspGrid = false;

template<int dim, class ctype>
inline constexpr bool isTensorProductYaspGrid<Dune::YaspGrid<dim, Dune::TensorProductCoordinates<ctype, dim>>> = true;

/*!
 * \brief The tissue grid, sequential for a self communicator and distributed otherwise
 *
 * The cells are equidistant, or graded symmetrically towards the boundary with a power law if the grid
 * has tensor-product coordinates. Both grids have faces in the plane x = 0.5 for an even number of cells.
 */
template<class Grid>
std::unique_ptr<Grid> makeBulkGrid(Dune::MPIHelper::MPICommunicator mpiComm)
{
    using Position = Dune::FieldVector<double, 3>;
    const auto lowerLeft = getParamFromGroup<Position>("Tissue", "Grid.LowerLeft");
    const auto upperRight = getParamFromGroup<Position>("Tissue", "Grid.UpperRight");
    const auto cells = getParamFromGroup<std::array<int, 3>>("Tissue", "Grid.Cells");
    const auto overlap = getParamFromGroup<int>("Tissue", "Grid.Overlap");
    const typename Grid::Communication comm(mpiComm);

    if constexpr (isTensorProductYaspGrid<Grid>)
    {
        const std::array<double, 3> exponents{1.3, 1.5, 0.8};
        std::array<std::vector<double>, 3> coordinates;
        for (int d = 0; d < 3; ++d)
        {
            const int half = cells[d]/2;
            coordinates[d].resize(cells[d] + 1);
            for (int i = 0; i <= half; ++i)
            {
                const auto t = 0.5*std::pow(double(i)/half, exponents[d]);
                coordinates[d][i] = lowerLeft[d] + t*(upperRight[d] - lowerLeft[d]);
                coordinates[d][cells[d] - i] = upperRight[d] - t*(upperRight[d] - lowerLeft[d]);
            }
        }
        return std::make_unique<Grid>(coordinates, std::bitset<3>(), overlap, comm);
    }
    else
        return std::make_unique<Grid>(lowerLeft, upperRight, cells, std::bitset<3>(), overlap, comm);
}


/*!
 * \brief Distribute the network over the processes of the tissue grid view, by the tissue partition or by
 *        a recursive coordinate bisection (parameter Network.Partition)
 */
template<class LowDimGrid, class BulkGridView>
void distributeNetwork(LowDimGrid& lowDimGrid, const BulkGridView& bulkGridView)
{
    lowDimGrid.setCommunicator(bulkGridView.comm());
    const auto partition = getParam<std::string>("Network.Partition");
    const auto owners = partition == "Tissue"
        ? lowDimPartition(bulkGridView, lowDimGrid.leafGridView())
        : Dune::FoamGridParallel::recursiveCoordinateBisection(lowDimGrid.leafGridView(), bulkGridView.comm().size());
    lowDimGrid.loadBalance(owners, getParamFromGroup<int>("Network", "Grid.Overlap"));
}

} // end namespace Dumux

#endif
