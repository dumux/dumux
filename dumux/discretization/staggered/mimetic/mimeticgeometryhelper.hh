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
 *        for the staggered discretization method
 */
#ifndef DUMUX_DISCRETIZATION_MIMETIC_GEOMETRY_HELPER_HH
#define DUMUX_DISCRETIZATION_MIMETIC_GEOMETRY_HELPER_HH

namespace Dumux
{

template<class GridView>
class MimeticGeometryHelper
{
    using Scalar = typename GridView::ctype;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using Intersection = typename GridView::Intersection;
    static constexpr int codimIntersection =  1;


public:
    MimeticGeometryHelper(const Intersection& intersection, const GridView& gridView)
    : intersection_(intersection), gridView_(gridView)
    {
    }

    /*!
    * \brief Returns the global dofIdx of the intersection itself
    */
   int dofIndex() const
   {
       //TODO: use proper intersection mapper!
       const auto inIdx = intersection_.indexInInside();
       return gridView_.indexSet().subIndex(intersection_.inside(), inIdx, codimIntersection);
   }

protected:

    // TODO: check whether to use references here or not
    const Intersection intersection_; //! The intersection of interest
    const GridView gridView_;
};




} // end namespace Dumux

#endif
