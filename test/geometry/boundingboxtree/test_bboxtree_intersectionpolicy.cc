//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Geometry
 * \brief Test that intersectingEntities uses the given intersection policy for all entities of the trees
 *
 * Vertical and horizontal segments cross in exactly one point each and do not overlap.
 * With the point policy, every crossing is an intersection with this point as its only corner.
 * With the segment policy, the default for two segments, there is no intersection.
 */
#include <config.h>

#include <array>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/geometryintersection.hh>
#include <dumux/geometry/intersectingentities.hh>

template<int dimworld>
void testPointPolicy()
{
    using namespace Dumux;
    using Point = Dune::FieldVector<double, dimworld>;
    using Segment = Dune::MultiLinearGeometry<double, 1, dimworld>;
    using SegmentSet = GeometriesEntitySet<Segment>;
    using PointPolicy = IntersectionPolicy::PointPolicy<double, dimworld>;

    // the segments lie in the plane z = 0.5 in 3d
    const auto point = [](double x, double y)
    {
        Point p(0.5);
        p[0] = x; p[1] = y;
        return p;
    };

    const std::vector<double> xs{0.1, 0.3, 0.45, 0.7, 0.9};
    const std::vector<double> ys{0.2, 0.35, 0.6, 0.8};
    std::vector<Segment> vertical, horizontal;
    for (const double x : xs)
        vertical.emplace_back(Dune::GeometryTypes::line, std::vector<Point>{point(x, 0.0), point(x, 1.0)});
    for (const double y : ys)
        horizontal.emplace_back(Dune::GeometryTypes::line, std::vector<Point>{point(0.0, y), point(1.0, y)});

    const BoundingBoxTree<SegmentSet> verticalTree(std::make_shared<SegmentSet>(vertical));
    const BoundingBoxTree<SegmentSet> horizontalTree(std::make_shared<SegmentSet>(horizontal));

    // every pair (i, j) of entities has to intersect once, in the point crossing(i, j)
    const auto checkCrossings = [&](const auto& intersections, std::size_t numI, std::size_t numJ,
                                    const auto& crossing, const std::string& name)
    {
        if (intersections.size() != numI*numJ)
            DUNE_THROW(Dune::Exception, name << ": found " << intersections.size()
                                             << " intersections, expected " << numI*numJ);

        std::vector<int> found(numI*numJ, 0);
        for (const auto& intersection : intersections)
        {
            const auto i = intersection.first();
            const auto j = intersection.second();
            if (i >= numI || j >= numJ)
                DUNE_THROW(Dune::Exception, name << ": invalid entity indices " << i << ", " << j);
            ++found[i*numJ + j];

            const auto& corners = intersection.corners();
            if (corners.size() != 1)
                DUNE_THROW(Dune::Exception, name << ": intersection of " << i << " and " << j << " has "
                                                 << corners.size() << " corners, expected a point");
            const Point expected = crossing(i, j);
            if ((corners[0] - expected).two_norm() > 1e-14)
                DUNE_THROW(Dune::Exception, name << ": intersection of " << i << " and " << j << " at "
                                                 << corners[0] << ", expected " << expected);
        }

        for (std::size_t k = 0; k < found.size(); ++k)
            if (found[k] != 1)
                DUNE_THROW(Dune::Exception, name << ": intersection of " << k/numJ << " and " << k%numJ
                                                 << " found " << found[k] << " times");

        std::cout << name << ": " << intersections.size() << " crossings found" << std::endl;
    };

    const std::string world = " (dimworld " + std::to_string(dimworld) + ")";

    checkCrossings(intersectingEntities(verticalTree, horizontalTree, PointPolicy{}),
                   xs.size(), ys.size(), [&](auto i, auto j){ return point(xs[i], ys[j]); },
                   "vertical tree, horizontal tree" + world);

    checkCrossings(intersectingEntities(horizontalTree, verticalTree, PointPolicy{}),
                   ys.size(), xs.size(), [&](auto i, auto j){ return point(xs[j], ys[i]); },
                   "horizontal tree, vertical tree" + world);

    for (std::size_t i = 0; i < horizontal.size(); ++i)
        checkCrossings(intersectingEntities(horizontal[i], verticalTree, PointPolicy{}),
                       1, xs.size(), [&](auto, auto j){ return point(xs[j], ys[i]); },
                       "horizontal segment " + std::to_string(i) + ", vertical tree" + world);

    if (const auto n = intersectingEntities(verticalTree, horizontalTree).size(); n != 0)
        DUNE_THROW(Dune::Exception, "Default policy for two segment trees" << world << ": found " << n
                                    << " intersections of crossing segments, expected none");

    if (const auto n = intersectingEntities(horizontal[0], verticalTree).size(); n != 0)
        DUNE_THROW(Dune::Exception, "Default policy for a segment and a segment tree" << world << ": found " << n
                                    << " intersections of crossing segments, expected none");
}

int main()
{
    testPointPolicy<2>();
    testPointPolicy<3>();
    return 0;
}
