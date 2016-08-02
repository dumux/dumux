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
    using CornerList = std::vector<GlobalPosition>;

    MpfaGeometryHelper(const typename Element::Geometry& elemGeom) : gt(elemGeom.type()) {}

    //! get sub control volume face geometries of an intersection for dim = 2
    static CornerList getScvfCorners(const typename Intersection::Geometry& geometry, unsigned int localIdx)
    {
        if (localIdx == 0)
            return CornerList({geometry.center(), geometry.corner(0)});
        else if (localIdx == 1)
            return CornerList({geometry.center(), geometry.corner(1)});
        else
            DUNE_THROW(Dune::InvalidStateException, "local index exceeds the number of corners of 2d intersections");
    }

    static GlobalPosition getScvfIntegrationPoint(const CornerList& scvfCorners, Scalar q)
    {
        auto d = scvfCorners[1];
        auto ip = scvfCorners[0];
        d -= ip;
        d *= q;
        ip += d;

        return ip;
    }

    static Scalar getScvfArea(const CornerList& scvfCorners)
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

} // end namespace

#endif