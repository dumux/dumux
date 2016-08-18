// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Base class for the finite volume geometry vector for box models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_GEOMETRYHELPER_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_GEOMETRYHELPER_HH

#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/math.hh>

namespace Dumux
{
//! A class to create sub control volume face geometries per intersection
template <class GridView, int dim>
class MpfaGeometryHelper
{};

//! Specialization for dim == 2
template <class GridView>
class MpfaGeometryHelper<GridView, 2>
{
private:
    using Scalar = typename GridView::ctype;
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using GlobalPosition = Dune::FieldVector<Scalar, dim>;
    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    using FaceReferenceElements = typename Dune::ReferenceElements<Scalar, dim-1>;

public:
    using PointVector = std::vector<GlobalPosition>;

    MpfaGeometryHelper(const typename Element::Geometry& elemGeom) : gt(elemGeom.type()) {}

    //! get sub control volume face corners of an intersection for the given local index
    static PointVector getScvfCorners(const typename Intersection::Geometry& geometry,
                                      unsigned int indexOnIntersection)
    {
        if (indexOnIntersection == 0)
            return PointVector({geometry.center(), geometry.corner(0)});
        else if (indexOnIntersection == 1)
            return PointVector({geometry.center(), geometry.corner(1)});
        else
            DUNE_THROW(Dune::InvalidStateException, "local index exceeds the number of corners of 2d intersections");
    }

    static GlobalPosition getScvfIntegrationPoint(const PointVector& scvfCorners, Scalar q)
    {
        auto d = scvfCorners[1];
        auto ip = scvfCorners[0];
        d -= ip;
        d *= q;
        ip += d;

        return ip;
    }

    static Scalar getScvfArea(const PointVector& scvfCorners)
    {
        return (scvfCorners[1]-scvfCorners[0]).two_norm();
    }

    std::size_t getNumLocalScvfs()
    {
        Dune::GeometryType triangle, quadrilateral;
        triangle.makeTriangle();
        quadrilateral.makeQuadrilateral();

        if (gt == triangle)
            return 6;
        else if (gt == quadrilateral)
            return 8;
        else
            DUNE_THROW(Dune::InvalidStateException, "unknown 2d geometry type " << gt);
    }

private:
    const Dune::GeometryType gt;
};

//! Specialization for dim == 3
template <class GridView>
class MpfaGeometryHelper<GridView, 3>
{
private:
    using Scalar = typename GridView::ctype;
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using GlobalPosition = Dune::FieldVector<Scalar, dim>;
    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    using FaceReferenceElements = typename Dune::ReferenceElements<Scalar, dim-1>;

public:
    using PointVector = std::vector<GlobalPosition>;

    MpfaGeometryHelper(const typename Element::Geometry& elemGeom) : gt(elemGeom.type()) {}

    //! get sub control volume face corners of an intersection for the given local index
    static PointVector getScvfCorners(const typename Intersection::Geometry& geometry,
                                      unsigned int indexOnIntersection)
    {
        // extract the corners of the sub control volumes
        const auto& referenceElement = FaceReferenceElements::general(geometry.type());

        // maximum number of necessary points is 9 (for quadrilateral)
        GlobalPosition p[9];
        auto corners = geometry.corners();

        // the intersection center
        p[0] = geometry.center();

        // vertices
        for (int i = 0; i < corners; ++i)
            p[i+1] = geometry.corner(i);

        // edge midpoints
        for (int i = 0; i < referenceElement.size(1); ++i)
            p[i+corners+1] = geometry.global(referenceElement.position(i, 1));

        // proceed according to number of corners
        switch (corners)
        {
        case 3: // triangle
        {
            //! Only build the maps the first time we encounter a triangle
            static const std::uint8_t vo = 1; //! vertex offset in point vector p
            static const std::uint8_t eo = 4; //! edge offset in point vector p
            static const std::uint8_t map[3][4] =
            {
                {0, eo+1, eo+0, vo+0},
                {0, eo+0, eo+2, vo+1},
                {0, eo+2, eo+1, vo+2}
            };

            return PointVector( {p[map[indexOnIntersection][0]],
                                 p[map[indexOnIntersection][1]],
                                 p[map[indexOnIntersection][2]],
                                 p[map[indexOnIntersection][3]]} );
        }
        case 4: // quadrilateral
        {
            //! Only build the maps the first time we encounter a quadrilateral
            static const std::uint8_t vo = 1; //! vertex offset in point vector p
            static const std::uint8_t eo = 5; //! face offset in point vector p
            static const std::uint8_t map[4][4] =
            {
                {0, eo+0, eo+2, vo+0},
                {0, eo+2, eo+1, vo+1},
                {0, eo+3, eo+0, vo+2},
                {0, eo+1, eo+3, vo+3}
            };

            return PointVector( {p[map[indexOnIntersection][0]],
                                 p[map[indexOnIntersection][1]],
                                 p[map[indexOnIntersection][2]],
                                 p[map[indexOnIntersection][3]]} );
        }
        default:
            DUNE_THROW(Dune::NotImplemented, "Box scvf boundary geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " corners=" << corners);
        }
    }

    static GlobalPosition getScvfIntegrationPoint(const PointVector& scvfCorners, Scalar q)
    {
        // in 3d, the integration point can not be moved from the midpoint
        return scvfCorners[0];
    }

    static Scalar getScvfArea(const PointVector& scvfCorners)
    {
        // after Wolfram alpha quadrilateral area
        return 0.5*Dumux::crossProduct(scvfCorners[3]-scvfCorners[0], scvfCorners[2]-scvfCorners[1]).two_norm();
    }

    std::size_t getNumLocalScvfs()
    {
        Dune::GeometryType tetrahedron, pyramid, prism, hexahedron;
        tetrahedron.makeTetrahedron();
        pyramid.makePyramid();
        prism.makePrism();
        hexahedron.makeHexahedron();

        if (gt == tetrahedron)
            return 12;
        else if (gt == pyramid)
            return 16;
        else if (gt == prism)
            return 18;
        else if (gt == hexahedron)
            return 24;
        else
            DUNE_THROW(Dune::InvalidStateException, "unknown 3d geometry type " << gt);
    }

private:
    const Dune::GeometryType gt; // the geometry type of the element
};

} // end namespace

#endif