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
#include <dumux/discretization/staggered/freeflow/staggeredgeometryhelper.hh>
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

          dofIdx_ = geometryHelper.dofIndex();
          oppositeIdx_ = geometryHelper.dofIndexOpposingFace();
          selfToOppositeDistance_ = geometryHelper.selfToOppositeDistance();

          pairData_ = geometryHelper.pairData();
          localFaceIdx_ = geometryHelper.localFaceIndex();
          dirIdx_ = geometryHelper.directionIndex();
          normalInPosCoordDir_ = unitOuterNormal()[directionIndex()] > 0.0;
          outerNormalScalar_ = unitOuterNormal()[directionIndex()];
          isGhostFace_ = false;
      }

      //! Constructor for a ghost face outside of the domain. Only needed to retrieve the center and scvIndices
      StaggeredSubControlVolumeFace(const GlobalPosition& dofPosition,
                                    const std::vector<IndexType>& scvIndices)
      {
          isGhostFace_ = true;
          center_ = dofPosition;
          scvIndices_ = scvIndices;
          scvfIndex_ = -1;
          dofIdx_ = -1;
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

    //! The center of the sub control volume face
    GlobalPosition dofPosition() const
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

    //! The global index of the dof living on this face
    IndexType dofIndex() const
    {
        return dofIdx_;
    }

    //! The global index of the dof living on the opposing face
    IndexType dofIndexOpposingFace() const
    {
        return oppositeIdx_;
    }

    //! The local index of this sub control volume face
    IndexType localFaceIdx() const
    {
        return localFaceIdx_;
    }

    //! Returns the dirction index of the facet (0 = x, 1 = y, 2 = z)
    int directionIndex() const
    {
        return dirIdx_;
    }

    //! The global index of this sub control volume face
    Scalar selfToOppositeDistance() const
    {
        return selfToOppositeDistance_;
    }

    //! The returns whether the unitNormal of the face point in positive coordinate direction
    bool normalInPosCoordDir() const
    {
        return normalInPosCoordDir_;
    }

    Scalar outerNormalScalar() const
    {
        return outerNormalScalar_;
    }


    auto pairData(const int idx) const
    {
        return pairData_[idx];
    }

    auto& pairData() const
    {
        return pairData_;
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

    int dofIdx_;
    int oppositeIdx_;
    Scalar selfToOppositeDistance_;
    std::vector<StaggeredSubFace> subfaces_;
    std::array<PairData<Scalar, GlobalPosition>, numPairs> pairData_;
    int localFaceIdx_;
    int dirIdx_;
    bool normalInPosCoordDir_;
    Scalar outerNormalScalar_;
    bool isGhostFace_;

};



} // end namespace

#endif
