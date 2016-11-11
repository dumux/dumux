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
 *        for the staggered discretizazion method
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_GEOMETRY_HELPER_HH
#define DUMUX_DISCRETIZATION_STAGGERED_GEOMETRY_HELPER_HH

#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/math.hh>

namespace Dumux
{

//! Create sub control volumes and sub control volume face geometries
template<class GridView, int dim>
class StaggeredGeometryHelper;


//! A class to create sub control volume and sub control volume face geometries per element
template <class GridView>
class StaggeredGeometryHelper<GridView, 2>
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

    StaggeredGeometryHelper(const Intersection& intersection, const GridView& gridView)
    : intersection_(intersection), elementGeometry_(intersection.inside().geometry()), gridView_(gridView)//, corners_(geometry.corners())
    {
    }

    int dofIdxSelf() const
    {
        //TODO: use proper intersection mapper!
        const auto inIdx = intersection_.indexInInside();
        const int numElements = gridView_.size(0);
        return gridView_.indexSet().subIndex(intersection_.inside(), inIdx, dim-1) + numElements;
    }

    int dofIdxOpposite() const
    {
        //TODO: use proper intersection mapper!
        const auto inIdx = intersection_.indexInInside();
        const int numElements = gridView_.size(0);
        const int localOppositeIdx = (inIdx % 2) ? (inIdx - 1) : (inIdx + 1);
        return gridView_.indexSet().subIndex(intersection_.inside(), localOppositeIdx, dim-1) + numElements;
    }

private:
    const Intersection& intersection_;
    const typename Element::Geometry& elementGeometry_; //! Reference to the element geometry
    const GridView gridView_;
//     GlobalPosition p[maxPoints]; // the points needed for construction of the geometries
//     std::size_t corners_; // number of element corners

};

} // end namespace Dumux

#endif
