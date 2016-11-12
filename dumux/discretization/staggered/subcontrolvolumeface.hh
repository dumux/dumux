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
#ifndef DUMUX_DISCRETIZATION_STAGGERED_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_STAGGERED_SUBCONTROLVOLUMEFACE_HH

#include <utility>
#include <dune/common/fvector.hh>
#include <dumux/discretization/subcontrolvolumefacebase.hh>
#include <dumux/discretization/staggered/staggeredgeometryhelper.hh>
#include <dumux/common/optional.hh>

#include <typeinfo>

namespace Dumux
{

template<class G, typename I>
class StaggeredSubFace
{
    using Geometry = G;
    using IndexType = I;
    using Scalar = typename Geometry::ctype;

private:
    std::vector<std::pair<int,int>> velocityDofIdxPair_;
    std::vector<Scalar> distance_;
    std::pair<int,int> elementPair_;
    int commonVertexIdx_;
};


/*!
 * \ingroup Discretization
 * \brief Class for a sub control volume face in the box method, i.e a part of the boundary
 *        of a sub control volume we compute fluxes on. We simply use the base class here.
 */
template<class G, typename I>
class StaggeredSubControlVolumeFace : public SubControlVolumeFaceBase<StaggeredSubControlVolumeFace<G, I>, G, I>
{
    using ParentType = SubControlVolumeFaceBase<StaggeredSubControlVolumeFace<G, I>, G, I>;
    using Geometry = G;
    using IndexType = I;

    using Scalar = typename Geometry::ctype;
    static const int dim = Geometry::mydimension;
    static const int dimworld = Geometry::coorddimension;

    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;
    using LocalPosition = Dune::FieldVector<Scalar, dim>;

    using StaggeredSubFace = Dumux::StaggeredSubFace<G,I>;

    static constexpr int numPairs = (dimworld == 2) ? 2 : 4;


public:
    // the default constructor
    StaggeredSubControlVolumeFace() = default;

    //! Constructor with intersection
    template <class Intersection, class GeometryHelper>
    StaggeredSubControlVolumeFace(const Intersection& is,
                               const typename Intersection::Geometry& isGeometry,
                               IndexType scvfIndex,
                               const std::vector<IndexType>& scvIndices,
                               const GeometryHelper& geometryHelper
                           )
    : ParentType(),
      geomType_(isGeometry.type()),
      area_(isGeometry.volume()),
      center_(isGeometry.center()),
      unitOuterNormal_(is.centerUnitOuterNormal()),
      scvfIndex_(scvfIndex),
      scvIndices_(scvIndices),
      boundary_(is.boundary())
      {
          corners_.resize(isGeometry.corners());
          for (int i = 0; i < isGeometry.corners(); ++i)
              corners_[i] = isGeometry.corner(i);

          selfIdx_ = geometryHelper.dofIdxSelf();
          oppositeIdx_ = geometryHelper.dofIdxOpposite();

          pairData_ = geometryHelper.pairData();
      }

    /*//! The copy constrcutor
    StaggeredSubControlVolumeFace(const StaggeredSubControlVolumeFace& other) = delete;

    //! The move constrcutor
    StaggeredSubControlVolumeFace(StaggeredSubControlVolumeFace&& other) = default;

    //! The copy assignment operator
    StaggeredSubControlVolumeFace& operator=(const StaggeredSubControlVolumeFace& other) = delete;

    //! The move assignment operator
    // StaggeredSubControlVolumeFace& operator=(StaggeredSubControlVolumeFace&& other)
    // {
    //     // We want to use the default copy/move assignment.
    //     // But since geometry is not copy assignable :( we
    //     // have to construct it again
    //     geometry_.release();
    //     geometry_.emplace(other.geometry_.value());
    //     unitOuterNormal_ = std::move(other.unitOuterNormal_);
    //     scvfIndex_ = std::move(other.scvfIndex_);
    //     scvIndices_ = std::move(other.scvIndices_);
    //     boundary_ = std::move(other.boundary_);
    //     return *this;
    // }*/

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
    IndexType outsideScvIdx() const
    {
        return scvIndices_[1];
    }

    //! The global index of this sub control volume face
    IndexType index() const
    {
        return scvfIndex_;
    }

    GlobalPosition corner(unsigned int localIdx) const
    {
        assert(localIdx < corners_.size() && "provided index exceeds the number of corners");
        return corners_[localIdx];
    }

    //! The geometry of the sub control volume face
    const Geometry geometry() const
    {
        return Geometry(geomType_, corners_);
    }

    //! The global index of this sub control volume face
    IndexType dofIndexSelf() const
    {
        return selfIdx_;
    }

    //! The global index of this sub control volume face
    IndexType dofIndexOpposite() const
    {
        return oppositeIdx_;
    }


    auto pairData(const int idx) const
    {
        return pairData_[idx];
    }

private:
    Dune::GeometryType geomType_;
    std::vector<GlobalPosition> corners_;
    Scalar area_;
    GlobalPosition center_;
    GlobalPosition unitOuterNormal_;
    IndexType scvfIndex_;
    std::vector<IndexType> scvIndices_;
    bool boundary_;

    int selfIdx_;
    int oppositeIdx_;
    std::vector<StaggeredSubFace> subfaces_;
    std::array<PairData<Scalar>, numPairs> pairData_;

};



} // end namespace

#endif
