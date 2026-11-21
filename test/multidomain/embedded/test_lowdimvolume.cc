// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test of the volume of a network of tubes in the bulk elements
 *
 * Straight segments inside an element, on a face and on an edge of the bulk grid have known
 * volumes \f$ \pi r^2 l \f$ that are split equally among the bulk elements containing the segment.
 */
#include <config.h>

#include <array>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/multidomain/embedded/lowdimvolume.hh>

namespace Dumux {

using Pos = Dune::FieldVector<double, 3>;

bool isClose(double a, double b)
{ return std::abs(a - b) <= 1e-13*std::max(std::abs(a), std::abs(b)); }

} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;
    initialize(argc, argv);

    using BulkGrid = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>>;
    BulkGrid bulkGrid(Pos(0.0), Pos({2.0, 2.0, 1.0}), std::array<int, 3>{2, 2, 1});

    // the segment centers identify the segments, the radius differs per segment
    struct Segment { Pos a, b; double radius; };
    const std::vector<Segment> segments = {
        {Pos({1.0, 0.5, 0.0}), Pos({1.0, 0.5, 1.0}), 0.1},  // on the face x = 1
        {Pos({1.0, 1.0, 0.0}), Pos({1.0, 1.0, 1.0}), 0.05}, // on the edge x = 1, y = 1
        {Pos({0.5, 0.5, 0.0}), Pos({0.5, 0.5, 1.0}), 0.2},  // inside the element (0, 0)
        {Pos({0.5, 1.5, 0.5}), Pos({1.5, 1.5, 0.5}), 0.1},  // across the face x = 1 inside the elements (0, 1) and (1, 1)
        {Pos({0.25, 0.25, 0.1}), Pos({1.75, 1.5, 0.9}), 0.1} // oblique through several elements
    };

    using LowDimGrid = Dune::FoamGrid<1, 3>;
    Dune::GridFactory<LowDimGrid> factory;
    for (const auto& s : segments)
    {
        factory.insertVertex(s.a);
        factory.insertVertex(s.b);
    }
    for (unsigned int i = 0; i < segments.size(); ++i)
        factory.insertElement(Dune::GeometryTypes::line, {2*i, 2*i + 1});
    const auto lowDimGrid = factory.createGrid();

    using BulkGridGeometry = CCTpfaFVGridGeometry<BulkGrid::LeafGridView>;
    using LowDimGridGeometry = CCTpfaFVGridGeometry<LowDimGrid::LeafGridView>;
    const BulkGridGeometry bulkGridGeometry(bulkGrid.leafGridView());
    const LowDimGridGeometry lowDimGridGeometry(lowDimGrid->leafGridView());

    // radius and length by element index
    std::vector<double> radii(segments.size()), lengths(segments.size());
    for (const auto& element : elements(lowDimGridGeometry.gridView()))
    {
        const auto center = element.geometry().center();
        for (const auto& s : segments)
            if ((center - 0.5*(s.a + s.b)).two_norm() < 1e-12)
            {
                const auto eIdx = lowDimGridGeometry.elementMapper().index(element);
                radii[eIdx] = s.radius;
                lengths[eIdx] = (s.b - s.a).two_norm();
            }
    }
    const auto radius = [&](std::size_t eIdx) { return radii[eIdx]; };

    const auto volume = EmbeddedCoupling::lowDimVolumePerBulkElement(bulkGridGeometry, lowDimGridGeometry, radius);

    int failures = 0;
    const auto check = [&](bool ok, const std::string& what)
    {
        if (!ok)
        {
            std::cout << "Failed: " << what << std::endl;
            ++failures;
        }
    };

    // the sum of the volumes is the sum of the volumes of the cylinders
    double expectedTotal = 0.0;
    for (const auto& s : segments)
        expectedTotal += M_PI*s.radius*s.radius*(s.b - s.a).two_norm();
    double total = 0.0;
    for (const auto v : volume)
        total += v;
    std::cout << "total volume " << total << ", sum of the cylinder volumes " << expectedTotal << std::endl;
    check(isClose(total, expectedTotal), "total volume");

    // element index i + 2*j of the cell (i, j)
    const auto cylinder = [&](int segment, double length) { return M_PI*segments[segment].radius*segments[segment].radius*length; };
    std::array<double, 4> expected{};
    expected[0] += cylinder(0, 1.0)/2.0 + cylinder(1, 1.0)/4.0 + cylinder(2, 1.0);
    expected[1] += cylinder(0, 1.0)/2.0 + cylinder(1, 1.0)/4.0;
    expected[2] += cylinder(1, 1.0)/4.0 + cylinder(3, 0.5);
    expected[3] += cylinder(1, 1.0)/4.0 + cylinder(3, 0.5);

    // the oblique segment: the part in each element is determined from the entry and exit points
    // of the straight line in the cells, which are computed independently of the glue
    {
        const auto& s = segments[4];
        const auto direction = s.b - s.a;
        const auto length = direction.two_norm();
        // the line crosses x = 1 at t = (1 - 0.25)/1.5 and y = 1 at t = (1 - 0.25)/1.25
        std::vector<double> ts = {0.0, 1.0, (1.0 - s.a[0])/direction[0], (1.0 - s.a[1])/direction[1]};
        std::sort(ts.begin(), ts.end());
        for (std::size_t i = 0; i + 1 < ts.size(); ++i)
        {
            const auto mid = s.a + 0.5*(ts[i] + ts[i+1])*direction;
            const int cell = (mid[0] > 1.0 ? 1 : 0) + 2*(mid[1] > 1.0 ? 1 : 0);
            expected[cell] += cylinder(4, (ts[i+1] - ts[i])*length);
        }
    }

    for (int i = 0; i < 4; ++i)
    {
        std::cout << "element " << i << ": " << volume[i] << ", expected " << expected[i] << std::endl;
        check(isClose(volume[i], expected[i]), "volume of element " + std::to_string(i));
    }

    // the glue-based overload gives the same result
    const auto glue = makeGlue(bulkGridGeometry, lowDimGridGeometry);
    const auto volumeFromGlue = EmbeddedCoupling::lowDimVolumePerBulkElement(glue, bulkGridGeometry, lowDimGridGeometry, radius);
    check(volumeFromGlue == volume, "glue overload");

    return failures == 0 ? 0 : 1;
}
