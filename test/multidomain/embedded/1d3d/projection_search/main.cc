// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test the closest segment surface search of the projection coupling manager
 *
 * The search through the bounding box tree has to return the same segment and the same
 * distance as a loop over all segments, including ties between segments at equal
 * signed distance (the segment with the lower index wins).
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <array>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d_projection.hh>

namespace {

using Grid = Dune::FoamGrid<1, 3>;
using GridView = Grid::LeafGridView;
using GridGeometry = Dumux::CCTpfaFVGridGeometry<GridView, true>;
using Network = Dumux::Detail::SegmentNetwork<GridGeometry>;
using Position = Dune::FieldVector<double, 3>;

struct Segments
{
    std::vector<std::array<Position, 2>> ends;
    std::vector<double> radii;
};

// random tree-like network with junctions, mirrored and duplicated segments
Segments makeSegments(std::size_t numSegments, bool variableRadius, std::mt19937_64& gen)
{
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::uniform_real_distribution<double> radius(0.002, 0.08);
    const double constantRadius = 0.03;

    std::vector<Position> nodes{ Position{ 0.3, 0.5, 0.5 } };
    Segments segments;
    auto addSegment = [&](const Position& a, const Position& b, double r)
    {
        segments.ends.push_back({ a, b });
        segments.radii.push_back(r);
    };

    while (segments.ends.size() < numSegments)
    {
        const double r = variableRadius ? radius(gen) : constantRadius;
        const auto from = nodes[std::uniform_int_distribution<std::size_t>(0, nodes.size()-1)(gen)];
        Position to = from;
        for (int i = 0; i < 3; ++i)
            to[i] += 0.3*(unit(gen) - 0.5);
        nodes.push_back(to);
        addSegment(from, to, r);

        const double kind = unit(gen);
        if (kind < 0.1 && segments.ends.size() < numSegments)
        {
            // mirror image on the plane x = 0 (query points on that plane are tied)
            auto mirror = [](Position p) { p[0] = -p[0]; return p; };
            addSegment(mirror(from), mirror(to), r);
        }
        else if (kind < 0.15 && segments.ends.size() < numSegments)
            addSegment(from, to, r);
    }
    return segments;
}

std::unique_ptr<Grid> makeGrid(const Segments& segments)
{
    Dune::GridFactory<Grid> factory;
    unsigned int vIdx = 0;
    for (const auto& [a, b] : segments.ends)
    {
        factory.insertVertex(a);
        factory.insertVertex(b);
        factory.insertElement(Dune::GeometryTypes::line, { vIdx, vIdx+1 });
        vIdx += 2;
    }
    return std::unique_ptr<Grid>(factory.createGrid());
}

// signed distance to the cylinder surface, evaluated independently of the network class
double signedDistance(const Position& p, const Position& a, const Position& b, double r)
{
    const auto ab = b - a;
    const auto ap = p - a;
    const double t = std::clamp((ab*ap)/(ab*ab), 0.0, 1.0);
    auto x = a;
    x.axpy(t, ab);
    return (p - x).two_norm() - r;
}

// loop over all segments, ties are resolved in favor of the lowest index
std::tuple<std::size_t, double> bruteForce(const Position& p, const Segments& segments)
{
    std::size_t best = 0;
    double bestDist = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < segments.ends.size(); ++i)
    {
        const double d = signedDistance(p, segments.ends[i][0], segments.ends[i][1], segments.radii[i]);
        if (d < bestDist)
        {
            bestDist = d;
            best = i;
        }
    }
    return { best, std::abs(bestDist) };
}

} // end namespace

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);

    std::mt19937_64 gen(20240611);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::size_t numFailures = 0;
    std::size_t numChecked = 0;
    std::size_t numTies = 0;

    auto fail = [&](const std::string& what)
    {
        if (++numFailures < 20)
            std::cerr << "FAILED: " << what << std::endl;
    };

    for (const std::size_t numSegments : { 1, 2, 7, 60, 800 })
    for (const bool variableRadius : { false, true })
    {
        const auto segments = makeSegments(numSegments, variableRadius, gen);
        const auto grid = makeGrid(segments);
        const GridGeometry gridGeometry(grid->leafGridView());
        const auto& mapper = gridGeometry.elementMapper();

        // the loop over all segments uses the element index of the grid
        Segments indexed = segments;
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            const auto eIdx = mapper.index(element);
            const auto g = element.geometry();
            indexed.ends[eIdx] = { g.corner(0), g.corner(1) };
            indexed.radii[eIdx] = segments.radii[eIdx];
        }

        const Network network(gridGeometry, [&](std::size_t eIdx) { return indexed.radii[eIdx]; });

        // query points: anywhere in and around the network, on the tube surfaces, inside the
        // tubes, far away, and on the plane of mirror symmetry
        std::vector<Position> points;
        for (int i = 0; i < 600; ++i)
        {
            Position p;
            const int kind = i % 6;
            for (int d = 0; d < 3; ++d)
                p[d] = unit(gen)*1.4 - 0.7 + 0.5*(d > 0);
            if (kind == 1)
                p[0] = 0.0;
            else if (kind == 2)
                p *= 50.0;
            else if (kind >= 3)
            {
                const auto s = std::uniform_int_distribution<std::size_t>(0, indexed.ends.size()-1)(gen);
                p = indexed.ends[s][0];
                p.axpy(unit(gen), indexed.ends[s][1] - indexed.ends[s][0]);
                const double offset = kind == 3 ? 0.0 : (kind == 4 ? 0.5 : 2.0)*indexed.radii[s]*unit(gen);
                p[1] += offset;
                p[2] += 1e-3*unit(gen);
            }
            points.push_back(p);
        }

        for (const auto& p : points)
        {
            const auto [bIdx, bDist] = bruteForce(p, indexed);
            const auto [tIdx, tDist] = network.findClosestSegmentSurface(p);
            ++numChecked;
            if (tIdx != bIdx || std::abs(tDist - bDist) > 1e-13)
                fail("N=" + std::to_string(numSegments) + " tree search (" + std::to_string(tIdx) + ", "
                     + std::to_string(tDist) + ") differs from loop (" + std::to_string(bIdx) + ", "
                     + std::to_string(bDist) + ")");

            // number of segments with the same signed distance as the closest one
            std::size_t numEqual = 0;
            for (std::size_t i = 0; i < indexed.ends.size(); ++i)
                numEqual += (signedDistance(p, indexed.ends[i][0], indexed.ends[i][1], indexed.radii[i])
                             == signedDistance(p, indexed.ends[bIdx][0], indexed.ends[bIdx][1], indexed.radii[bIdx]));
            numTies += (numEqual > 1);
        }

        // the candidate list is sorted, valid for all points within the margin, and selective
        std::size_t numCandidates = 0;
        const double margin = 0.05;
        for (const auto& p : points)
        {
            const auto candidates = network.closestSegmentSurfaceCandidates(p, margin);
            numCandidates += candidates.size();
            if (!std::is_sorted(candidates.begin(), candidates.end()))
                fail("candidates not sorted");

            for (int k = 0; k < 8; ++k)
            {
                Position q = p;
                Position dir;
                for (int d = 0; d < 3; ++d)
                    dir[d] = unit(gen) - 0.5;
                dir *= margin*(k == 7 ? 1.0 : unit(gen))/dir.two_norm();
                q += dir;

                const auto [bIdx, bDist] = bruteForce(q, indexed);
                const auto [tIdx, tDist] = network.findClosestSegmentSurface(q, candidates);
                ++numChecked;
                if (tIdx != bIdx || std::abs(tDist - bDist) > 1e-13)
                    fail("N=" + std::to_string(numSegments) + " search among face candidates differs from loop");
            }
        }

        std::cout << "N=" << numSegments << (variableRadius ? " variable" : " constant")
                  << " radius: mean candidate count with margin " << double(numCandidates)/points.size()
                  << " of " << numSegments << std::endl;

        if (numSegments >= 800 && double(numCandidates)/points.size() > 0.25*numSegments)
            fail("candidate list is not selective");
    }

    if (numTies == 0)
        fail("no query point with tied segments was generated");

    std::cout << numChecked << " searches, " << numTies << " with ties, " << numFailures << " failures" << std::endl;
    return numFailures == 0 ? 0 : 1;
}
