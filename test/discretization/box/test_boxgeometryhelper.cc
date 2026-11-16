// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test the static box geometry helper functions, which construct the sub-control-volume
 *        (face) corners from a geometry type and a transformation, against the functions
 *        bound to an element geometry
 *
 * The corners obtained with the transformation to global coordinates have to coincide with
 * the corners of the bound helper, and the local scvf centers mapped to global coordinates
 * with the centers of the scvf corners.
 */
#include <config.h>

#include <iostream>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/geometry/center.hh>
#include <dumux/discretization/box/boxgeometryhelper.hh>
#include <dumux/discretization/pq1bubble/geometryhelper.hh>

namespace Mock {

template<int dim, int dimWorld>
struct GridView
{
    template<int i>
    struct Codim
    {
        struct Entity { using Geometry = Dune::MultiLinearGeometry<double, dim, dimWorld>; };
    };

    struct IndexSet
    {
        using IndexType = int;
    };

    struct Intersection { using Geometry = int; };

    using ctype = double;
    static constexpr int dimension = dim;
    static constexpr int dimensionworld = dimWorld;
};

template<int dim, int dimWorld>
struct SubGeoType
{
    struct Traits
    {
        using CornerStorage = std::vector<Dune::FieldVector<double, dimWorld>>;
        using Geometry = Dune::MultiLinearGeometry<double, dim, dimWorld>;
        using LocalIndexType = int;
    };
};

} // end namespace Mock

namespace {

template<class Corners>
void checkCorners(const Corners& expected, const Corners& computed, const std::string& what)
{
    if (expected.size() != computed.size())
        DUNE_THROW(Dune::Exception, what << ": " << computed.size() << " corners instead of " << expected.size());

    for (std::size_t i = 0; i < expected.size(); ++i)
        if ((expected[i] - computed[i]).two_norm() > 1e-13)
            DUNE_THROW(Dune::Exception, what << ": corner " << i << " is " << computed[i] << " instead of " << expected[i]);
}

template<class Position>
void checkPosition(const Position& expected, const Position& computed, const std::string& what)
{
    if ((expected - computed).two_norm() > 1e-13)
        DUNE_THROW(Dune::Exception, what << ": " << computed << " instead of " << expected);
}

// The element is the image of the reference element under an affine map,
// so that local scvf centers map to the centers of the global scvf corners.
template<int dim, int dimWorld>
Dune::MultiLinearGeometry<double, dim, dimWorld> makeAffineGeometry(Dune::GeometryType type)
{
    const auto& ref = Dune::referenceElement<double, dim>(type);
    std::vector<Dune::FieldVector<double, dimWorld>> corners;
    for (int i = 0; i < ref.size(dim); ++i)
    {
        const auto local = ref.position(i, dim);
        Dune::FieldVector<double, dimWorld> corner;
        for (int k = 0; k < dimWorld; ++k)
        {
            corner[k] = 0.1*(k + 1);
            for (int l = 0; l < dim; ++l)
                corner[k] += (k == l ? 1.0 + 0.5*k : 0.3 - 0.1*l)*local[l];
        }
        corners.push_back(corner);
    }

    return { type, corners };
}

template<int dim, int dimWorld>
void testBoxGeometryHelper(Dune::GeometryType type)
{
    std::cout << "Testing " << type << " in " << dimWorld << "d" << std::endl;

    using GridView = Mock::GridView<dim, dimWorld>;
    using ScvType = Mock::SubGeoType<dim, dimWorld>;
    using ScvfType = Mock::SubGeoType<dim-1, dimWorld>;
    using Helper = Dumux::BoxGeometryHelper<GridView, dim, ScvType, ScvfType>;

    const auto geometry = makeAffineGeometry<dim, dimWorld>(type);
    const auto& ref = Dune::referenceElement<double, dim>(type);
    const auto toGlobal = [&](const auto& local){ return geometry.global(local); };
    const Helper helper(geometry);

    for (int i = 0; i < ref.size(dim); ++i)
        checkCorners(helper.getScvCorners(i), Helper::getScvCorners(type, toGlobal, i), "scv corners");

    for (unsigned int i = 0; i < helper.numInteriorScvf(); ++i)
    {
        const auto scvfCorners = helper.getScvfCorners(i);
        checkCorners(scvfCorners, Helper::getScvfCorners(type, toGlobal, i), "scvf corners");

        if constexpr (dim == dimWorld)
            checkPosition(Dumux::center(scvfCorners), geometry.global(Helper::localScvfCenter(type, i)), "scvf center");
    }

    for (int facet = 0; facet < ref.size(1); ++facet)
    {
        for (int i = 0; i < ref.size(facet, 1, dim); ++i)
        {
            const auto scvfCorners = helper.getBoundaryScvfCorners(facet, i);
            checkCorners(scvfCorners, Helper::getBoundaryScvfCorners(type, toGlobal, facet, i), "boundary scvf corners");

            if constexpr (dim == dimWorld)
            {
                using PQ1BubbleHelper = Dumux::PQ1BubbleGeometryHelper<GridView, ScvType, ScvfType>;
                const auto center = Dumux::center(scvfCorners);
                checkPosition(center, geometry.global(Helper::localBoundaryScvfCenter(type, facet, i)), "boundary scvf center");
                checkPosition(center, geometry.global(PQ1BubbleHelper::localBoundaryScvfCenter(type, facet, i)), "pq1bubble boundary scvf center");
            }
        }
    }
}

} // end namespace

int main()
{
    using namespace Dune::GeometryTypes;

    testBoxGeometryHelper<1, 1>(line);
    testBoxGeometryHelper<1, 2>(line);
    testBoxGeometryHelper<1, 3>(line);
    testBoxGeometryHelper<2, 2>(triangle);
    testBoxGeometryHelper<2, 2>(quadrilateral);
    testBoxGeometryHelper<2, 3>(triangle);
    testBoxGeometryHelper<2, 3>(quadrilateral);
    testBoxGeometryHelper<3, 3>(tetrahedron);
    testBoxGeometryHelper<3, 3>(hexahedron);
    testBoxGeometryHelper<3, 3>(prism);

    return 0;
}
