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
#ifndef DUMUX_DISCRETIZATION_MIMETIC_CP_GEOMETRY_HELPER_HH
#define DUMUX_DISCRETIZATION_MIMETIC_CP_GEOMETRY_HELPER_HH

#include <dumux/common/math.hh>

namespace Dumux
{

template<class GridView>
class MimeticCPGeometryHelper
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
    MimeticCPGeometryHelper(const Element& element, const GridView& gridView)
        : elementGeometry_(element.geometry()),
          gridView_(gridView),
          dofIndex_(0),
          localIndex_(-1),
          area_(0),
          unitOuterNormal_(0)
    {
        calculateIntegratedNormals();
        totalFaceArea_.resize(6,0.0);

        // fill neighbor information and control volume face data:
        for (const auto& intersection : Dune::intersections(gridView, element))
        {
            totalFaceArea_[intersection.indexInInside()] += intersection.geometry().volume();
        }
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
  int localFaceIndex() const
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

  template<class IntersectionMapper>
  void updateLocalFace(const IntersectionMapper& intersectionMapper_, const Intersection& intersection)
  {
      localIndex_++;
      dofIndex_ = intersectionMapper_.globalIntersectionIndex(intersection.inside(), localIndex_);;
      auto integratedNormal = integratedNormals_[intersection.indexInInside()];
      Scalar integratedNormalNorm = integratedNormal.two_norm();
      unitOuterNormal_ = integratedNormal;
      unitOuterNormal_ /= integratedNormalNorm;
      area_ = integratedNormalNorm;
      //For non-matching grids multiplication with ratio
      area_ *= intersection.geometry().volume()/totalFaceArea_[intersection.indexInInside()];
  }

protected:

  void calculateIntegratedNormals()
  {
      if(elementGeometry_.corners() != 8)
          DUNE_THROW(Dune::NotImplemented, "Integration of Normals only implemented for hexahedral grids!");

      unsigned int numFaces = 6;
      integratedNormals_.resize(numFaces);

      for(int faceIdx=0; faceIdx < numFaces; faceIdx++)
      {
          GlobalPosition x1,x2,x3,x4;
          switch(faceIdx)
          {
              case 0:
                  x1 = elementGeometry_.corner(2);
                  x2 = elementGeometry_.corner(0);
                  x3 = elementGeometry_.corner(6);
                  x4 = elementGeometry_.corner(4);
                  break;
              case 1:
                  x1 = elementGeometry_.corner(1);
                  x2 = elementGeometry_.corner(3);
                  x3 = elementGeometry_.corner(5);
                  x4 = elementGeometry_.corner(7);
                  break;
              case 2:
                  x1 = elementGeometry_.corner(0);
                  x2 = elementGeometry_.corner(1);
                  x3 = elementGeometry_.corner(4);
                  x4 = elementGeometry_.corner(5);
                  break;
              case 3:
                  x1 = elementGeometry_.corner(3);
                  x2 = elementGeometry_.corner(2);
                  x3 = elementGeometry_.corner(7);
                  x4 = elementGeometry_.corner(6);
                  break;
              case 4:
                  x1 = elementGeometry_.corner(2);
                  x2 = elementGeometry_.corner(3);
                  x3 = elementGeometry_.corner(0);
                  x4 = elementGeometry_.corner(1);
                  break;
              case 5:
                  x1 = elementGeometry_.corner(4);
                  x2 = elementGeometry_.corner(5);
                  x3 = elementGeometry_.corner(6);
                  x4 = elementGeometry_.corner(7);
                  break;
          }

          integratedNormals_[faceIdx] =  crossProduct(x1,x2) - crossProduct(x1,x3) + crossProduct(x2,x4) - crossProduct(x3,x4);
          integratedNormals_[faceIdx] *= 0.5;

      }
  }

    // TODO: check whether to use references here or not
    int dofIndex_; //! The global dofIdx of the intersectio
    int localIndex_; //! The local index of the intersection
    Scalar area_;
    GlobalPosition unitOuterNormal_;
    std::vector<GlobalPosition> integratedNormals_;
    std::vector<Scalar> totalFaceArea_;
    const typename Element::Geometry& elementGeometry_; //! Reference to the element geometry
    const GridView gridView_;
};




} // end namespace Dumux

#endif
