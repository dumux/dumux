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

    using Element = typename GridView::template Codim<0>::Entity;
    using Geometry = typename Element::Geometry;
    using Intersection = typename GridView::Intersection;
    static constexpr int codimIntersection =  1;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;


public:
    MimeticGeometryHelper(const Element& element, const GridView& gridView)
    :gridView_(gridView),
     dofIndex_(0),
     localIndex_(0),
     area_(0),
     unitOuterNormal_(0)
    {
    }

    /*!
    * \brief Returns the global dofIdx of the intersection itself
    */
   int dofIndex() const
   {
       return dofIndex_;
   }

   /*!
   * \brief Returns the local index of the intersection itself
   */
  int localIndex() const
  {
      return localIndex_;
  }

  //! The area of the sub control volume face
  Scalar area() const
  {
      return area_;
  }

  GlobalPosition unitOuterNormal() const
  {
      return unitOuterNormal_;
  }

  void updateLocalFace(const int dofIndex, const int localIndex, const Intersection& intersection)
  {
      dofIndex_ = dofIndex;
      localIndex_ = localIndex;
      unitOuterNormal_ = intersection.centerUnitOuterNormal();
      area_ = intersection.geometry().volume();
  }

protected:

    // TODO: check whether to use references here or not
    int dofIndex_; //! The global dofIdx of the intersectio
    int localIndex_; //! The local index of the intersection
    Scalar area_;
    GlobalPosition unitOuterNormal_;
    const GridView gridView_;
};




} // end namespace Dumux

#endif
