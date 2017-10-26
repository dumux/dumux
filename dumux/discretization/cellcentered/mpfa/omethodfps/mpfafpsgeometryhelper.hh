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

#include <dumux/common/math.hh>
#include <dumux/common/optional.hh>

namespace Dumux
{

//! A class to mainly extract the scv corners of the scvs inside an element
template<class GridView, int dim>
class MpfaFpsGeometryHelper;

//! Specialization for dim = 2
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

} // end namespace Dumux

#endif
