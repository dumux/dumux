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
 * \brief Base class for a sub control volume face
 */
#ifndef DUMUX_DISCRETIZATION_MIMETIC_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_MIMETIC_SUBCONTROLVOLUMEFACE_HH

#include <utility>
#include <dune/common/fvector.hh>
#include <dumux/discretization/subcontrolvolumefacebase.hh>
#include <dumux/common/optional.hh>

#include <typeinfo>
#include "mimeticgeometryhelper.hh"

namespace Dumux
{

/*!
 * \ingroup Discretization
 * \brief Class for a sub control volume face in the box method, i.e a part of the boundary
 *        of a sub control volume we compute fluxes on. We simply use the base class here.
 */
template<class G, typename I>
class MimeticSubControlVolumeFace : public SubControlVolumeFaceBase<MimeticSubControlVolumeFace<G, I>, G, I>
{
    using ParentType = SubControlVolumeFaceBase<MimeticSubControlVolumeFace<G, I>, G, I>;
    using Geometry = G;
    using IndexType = I;

    using Scalar = typename Geometry::ctype;
    static const int dim = Geometry::mydimension;
    static const int dimworld = Geometry::coorddimension;

    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;
    using LocalPosition = Dune::FieldVector<Scalar, dim>;

    static constexpr int numPairs = (dimworld == 2) ? 2 : 4;


public:
    // the default constructor
    MimeticSubControlVolumeFace() = default;

    //! Constructor with intersection
    template <class Intersection, class GeometryHelper>
    MimeticSubControlVolumeFace(const Intersection& is,
                               const typename Intersection::Geometry& isGeometry,
                               IndexType scvfIndex,
                               const std::vector<IndexType>& scvIndices,
                               const GeometryHelper& geometryHelper
                           )
    : ParentType(),
      area_(geometryHelper.area()),
      center_(isGeometry.center()),
      unitOuterNormal_(geometryHelper.unitOuterNormal()),
      scvfIndex_(scvfIndex),
      scvIndices_(scvIndices),
      boundary_(is.boundary())
      {
          dofIdx_ = geometryHelper.dofIndex();
          localFaceIdx_ = geometryHelper.localFaceIndex();
      }


    //! The center of the sub control volume face
    GlobalPosition center() const
    {
        return center_;
    }

    //! The integration point for flux evaluations in global coordinates
    GlobalPosition ipGlobal() const
    {
        // Return center for now
        return center_;
    }

    //! The area of the sub control volume face
    Scalar area() const
    {
        return area_;
    }

    //! returns bolean if the sub control volume face is on the boundary
    bool boundary() const
    {
        return boundary_;
    }

    GlobalPosition unitOuterNormal() const
    {
        return unitOuterNormal_;
    }

    //! index of the inside sub control volume for spatial param evaluation
    IndexType insideScvIdx() const
    {
        return scvIndices_[0];
    }

    //! index of the outside sub control volume for spatial param evaluation
    // This results in undefined behaviour if boundary is true
    IndexType outsideScvIdx(int i = 0) const
    {
        return scvIndices_[i+1];
    }

    //! The number of outside scvs connection via this scv face
    std::size_t numOutsideScvs() const
    {
        return scvIndices_.size()-1;
    }

    //! The global index of this sub control volume face
    IndexType index() const
    {
        return scvfIndex_;
    }

    //! The global index of the dof living on this face
    IndexType dofIndex() const
    {
        return dofIdx_;
    }

    //! The local index of this sub control volume face
    IndexType localFaceIdx() const
    {
        return localFaceIdx_;
    }

    int directionIndex() const
    {
        return 0;
    }

    int fluxMultiplier() const
    {
        if(scvIndices_[0] < scvIndices_[1])
            return 1;
        else
            return -1;
    }

private:
    Scalar area_;
    GlobalPosition center_;
    GlobalPosition unitOuterNormal_;
    IndexType scvfIndex_;
    std::vector<IndexType> scvIndices_;
    bool boundary_;
    int dofIdx_;
    int localFaceIdx_;
};



} // end namespace

#endif
