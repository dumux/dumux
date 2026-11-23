//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test that a bounding box tree of an empty entity set can be built and queried
 *
 * A process of a distributed grid may hold no elements of a grid, e.g. of a network that does not
 * reach its part of the domain. All queries of an empty tree find nothing.
 */
#include <config.h>

#include <array>
#include <limits>
#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/distance.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/intersectingentities.hh>

int main()
{
    using namespace Dumux;
    using Point = Dune::FieldVector<double, 3>;
    using Segment = Dune::MultiLinearGeometry<double, 1, 3>;
    using SegmentSet = GeometriesEntitySet<Segment>;

    const BoundingBoxTree<SegmentSet> emptyTree(std::make_shared<SegmentSet>(std::vector<Segment>{}));
    const Segment segment(Dune::GeometryTypes::line, std::vector<Point>{Point(0.0), Point(1.0)});
    const BoundingBoxTree<SegmentSet> tree(std::make_shared<SegmentSet>(std::vector<Segment>{segment}));

    if (emptyTree.numBoundingBoxes() != 0)
        DUNE_THROW(Dune::Exception, "The tree of an empty set has " << emptyTree.numBoundingBoxes() << " bounding boxes");

    if (!intersectingEntities(Point(0.5), emptyTree).empty())
        DUNE_THROW(Dune::Exception, "A point intersects an empty tree");
    if (!intersectingEntities(segment, emptyTree).empty())
        DUNE_THROW(Dune::Exception, "A segment intersects an empty tree");
    if (!intersectingEntities(emptyTree, tree).empty() || !intersectingEntities(tree, emptyTree).empty())
        DUNE_THROW(Dune::Exception, "A tree intersects an empty tree");

    const std::array<double, 6> box{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    if (!intersectingBoxes(box.data(), emptyTree).empty())
        DUNE_THROW(Dune::Exception, "A box intersects an empty tree");
    if (!intersectingBoxes(emptyTree, tree).empty() || !intersectingBoxes(tree, emptyTree).empty())
        DUNE_THROW(Dune::Exception, "The boxes of a tree intersect an empty tree");

    const auto [squaredDistance, index] = closestEntity(Point(0.5), emptyTree);
    if (squaredDistance != std::numeric_limits<double>::max())
        DUNE_THROW(Dune::Exception, "An empty tree has an entity at squared distance " << squaredDistance);

    return 0;
}
