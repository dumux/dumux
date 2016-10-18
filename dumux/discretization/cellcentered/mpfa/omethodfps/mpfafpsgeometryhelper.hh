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
 * \brief Helper class constructing the dual grid finite volume geometries
 *        for the box discretizazion method
 */
#ifndef DUMUX_DISCRETIZATION_MPFA_FPS_GEOMETRY_HELPER_HH
#define DUMUX_DISCRETIZATION_MPFA_FPS_GEOMETRY_HELPER_HH

#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dumux/common/math.hh>
#include <dumux/common/optional.hh>

namespace Dumux
{

//! Create sub control volumes and sub control volume face geometries
template<class GridView, int dim>
class MpfaFpsGeometryHelper;

//! A class to create sub control volume and sub control volume face geometries per element
template <class GridView>
class MpfaFpsGeometryHelper<GridView, 2>
{
    using Scalar = typename GridView::ctype;
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using ScvGeometry = Dune::CachedMultiLinearGeometry<Scalar, dim, dimWorld>;
    using ScvfGeometry = Dune::CachedMultiLinearGeometry<Scalar, dim-1, dimWorld>;

    using GlobalPosition = typename ScvGeometry::GlobalCoordinate;
    using PointVector = std::vector<GlobalPosition>;

    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    using ReferenceElements = typename Dune::ReferenceElements<Scalar, dim>;

    //! the maximum number of helper points used to construct the geometries
    //! Using a statically sized point array is much faster than dynamic allocation
    static constexpr int maxPoints = 9;

public:

    MpfaFpsGeometryHelper(const typename Element::Geometry& geometry)
    : corners_(geometry.corners())
    {
        // extract the corners of the sub control volumes
        const auto& referenceElement = ReferenceElements::general(geometry.type());

        // the element center
        p[0] = geometry.center();

        // vertices
        for (int i = 0; i < corners_; ++i)
            p[i+1] = geometry.corner(i);

        // face midpoints
        for (int i = 0; i < referenceElement.size(1); ++i)
            p[i+corners_+1] = geometry.global(referenceElement.position(i, 1));
    }

    //! Create a vector with the scv corners
    PointVector getScvCorners(unsigned int localScvIdx) const
    {
        // proceed according to number of corners of the element
        switch (corners_)
        {
        case 3: // triangle
        {
            //! Only build the maps the first time we encounter a triangle
            static const std::uint8_t vo = 1; //! vertex offset in point vector p
            static const std::uint8_t fo = 4; //! face offset in point vector p
            static const std::uint8_t map[3][4] =
            {
                {0, fo+1, fo+0, vo+0},
                {0, fo+0, fo+2, vo+1},
                {0, fo+2, fo+1, vo+2}
            };

            return PointVector( {p[map[localScvIdx][0]],
                                 p[map[localScvIdx][1]],
                                 p[map[localScvIdx][2]],
                                 p[map[localScvIdx][3]]} );
        }
        case 4: // quadrilateral
        {
            //! Only build the maps the first time we encounter a quadrilateral
            static const std::uint8_t vo = 1; //! vertex offset in point vector p
            static const std::uint8_t fo = 5; //! face offset in point vector p
            static const std::uint8_t map[4][4] =
            {
                {0, fo+0, fo+2, vo+0},
                {0, fo+2, fo+1, vo+1},
                {0, fo+3, fo+0, vo+2},
                {0, fo+1, fo+3, vo+3}
            };

            return PointVector( {p[map[localScvIdx][0]],
                                 p[map[localScvIdx][1]],
                                 p[map[localScvIdx][2]],
                                 p[map[localScvIdx][3]]} );
        }
        default:
            DUNE_THROW(Dune::NotImplemented, "Mpfa-Fps scv corners for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " corners=" << corners_);
        }
    }

private:
    GlobalPosition p[maxPoints]; // the points needed for construction of the geometries
    std::size_t corners_; // number of element corners
};

//! A class to create sub control volume and sub control volume face geometries per element
template <class GridView>
class MpfaFpsGeometryHelper<GridView, 3>
{
    using Scalar = typename GridView::ctype;
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using ScvGeometry = Dune::CachedMultiLinearGeometry<Scalar, dim, dimWorld>;
    using ScvfGeometry = Dune::CachedMultiLinearGeometry<Scalar, dim-1, dimWorld>;

    using GlobalPosition = typename ScvGeometry::GlobalCoordinate;
    using PointVector = std::vector<GlobalPosition>;

    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    using ReferenceElements = typename Dune::ReferenceElements<Scalar, dim>;
    using FaceReferenceElements = typename Dune::ReferenceElements<Scalar, dim-1>;

    //! the maximum number of helper points used to construct the geometries
    //! Using a statically sized point array is much faster than dynamic allocation
    static constexpr int maxPoints = 27;
public:
    MpfaFpsGeometryHelper(const typename Element::Geometry& geometry)
    : elementGeometry_(geometry), corners_(geometry.corners())
    {
        // extract the corners of the sub control volumes
        const auto& referenceElement = ReferenceElements::general(geometry.type());

        // the element center
        p[0] = geometry.center();

        // vertices
        for (int i = 0; i < corners_; ++i)
            p[i+1] = geometry.corner(i);

        // edge midpoints
        for (int i = 0; i < referenceElement.size(dim-1); ++i)
            p[i+corners_+1] = geometry.global(referenceElement.position(i, dim-1));

        // face midpoints
        for (int i = 0; i < referenceElement.size(1); ++i)
            p[i+corners_+1+referenceElement.size(dim-1)] = geometry.global(referenceElement.position(i, 1));
    }

    //! Create a vector with the scv corners
    PointVector getScvCorners(unsigned int localScvIdx) const
    {
        // proceed according to number of corners of the element
        switch (corners_)
        {
        case 4: // tetrahedron
        {
            //! Only build the maps the first time we encounter a tetrahedron
            static const std::uint8_t vo = 1; //! vertex offset in point vector p
            static const std::uint8_t eo = 5; //! edge offset in point vector p
            static const std::uint8_t fo = 11; //! face offset in point vector p
            static const std::uint8_t map[4][8] =
            {
                {vo+0, eo+0, eo+1, fo+0, eo+3, fo+1, fo+2,    0},
                {vo+1, eo+2, eo+0, fo+0, eo+4, fo+3, fo+1,    0},
                {vo+2, eo+1, eo+2, fo+0, eo+5, fo+2, fo+3,    0},
                {vo+3, eo+3, eo+5, fo+2, eo+4, fo+1, fo+3,    0}
            };

            return PointVector( {p[map[localScvIdx][0]],
                                 p[map[localScvIdx][1]],
                                 p[map[localScvIdx][2]],
                                 p[map[localScvIdx][3]],
                                 p[map[localScvIdx][4]],
                                 p[map[localScvIdx][5]],
                                 p[map[localScvIdx][6]],
                                 p[map[localScvIdx][7]]} );
        }
        case 8: // hexahedron
        {
            //! Only build the maps the first time we encounter a quadrilateral
            static const std::uint8_t vo = 1; //! vertex offset in point vector p
            static const std::uint8_t eo = 9; //! edge offset in point vector p
            static const std::uint8_t fo = 21; //! face offset in point vector p
            static const std::uint8_t map[8][8] =
            {
                {vo+0, eo+6, eo+4, fo+4, eo+0, fo+2, fo+0,    0},
                {vo+1, eo+5, eo+6, fo+4, eo+1, fo+1, fo+2,    0},
                {vo+2, eo+4, eo+7, fo+4, eo+2, fo+0, fo+3,    0},
                {vo+3, eo+7, eo+5, fo+4, eo+3, fo+3, fo+1,    0},
                {vo+4, eo+8, eo+10, fo+5, eo+0, fo+0, fo+2,    0},
                {vo+5, eo+10, eo+9, fo+5, eo+1, fo+2, fo+1,    0},
                {vo+6, eo+11, eo+8, fo+5, eo+2, fo+3, fo+0,    0},
                {vo+7, eo+9, eo+11, fo+5, eo+3, fo+1, fo+3,    0},
            };

            return PointVector( {p[map[localScvIdx][0]],
                                 p[map[localScvIdx][1]],
                                 p[map[localScvIdx][2]],
                                 p[map[localScvIdx][3]],
                                 p[map[localScvIdx][4]],
                                 p[map[localScvIdx][5]],
                                 p[map[localScvIdx][6]],
                                 p[map[localScvIdx][7]]} );
        }
        default:
            DUNE_THROW(Dune::NotImplemented, "Mpfa-Fps scv corners for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " corners=" << corners_);
        }
    }

private:
    const typename Element::Geometry& elementGeometry_; //! Reference to the element geometry
    GlobalPosition p[maxPoints]; // the points needed for construction of the scv/scvf geometries
    std::size_t corners_; // number of element corners
};

} // end namespace Dumux

#endif
