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
 * \brief Specilialized sub control volume face for free-flow staggered grid models
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FREE_FLOW_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_STAGGERED_FREE_FLOW_SUBCONTROLVOLUMEFACE_HH

#include <utility>
#include <dune/common/fvector.hh>

#include <dumux/discretization/subcontrolvolumefacebase.hh>
#include <dumux/discretization/staggered/freeflow/staggeredgeometryhelper.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/optional.hh>

#include <typeinfo>

namespace Dumux
{

/*!
 * \ingroup Discretization
 * \brief Class for a sub control volume face in the staggered method, i.e a part of the boundary
 *        of a sub control volume we compute fluxes on. This is a specialization for free flow models.
 */
template<class ScvfGeometryTraits>
class FreeFlowStaggeredSubControlVolumeFace : public SubControlVolumeFaceBase<FreeFlowStaggeredSubControlVolumeFace<ScvfGeometryTraits>, ScvfGeometryTraits>
{
    using ParentType = SubControlVolumeFaceBase<FreeFlowStaggeredSubControlVolumeFace<ScvfGeometryTraits>,ScvfGeometryTraits>;
    using Geometry = typename ScvfGeometryTraits::Geometry;
    using GridIndexType = typename ScvfGeometryTraits::GridIndexType;

    using Scalar = typename ScvfGeometryTraits::Scalar;
    static const int dim = Geometry::mydimension;
    static const int dimworld = Geometry::coorddimension;

    using GlobalPosition = typename ScvfGeometryTraits::GlobalPosition;

    static constexpr int numPairs = (dimworld == 2) ? 2 : 4;

public:
    //! state the traits public and thus export all types
    using Traits = ScvfGeometryTraits;

    // the default constructor
    FreeFlowStaggeredSubControlVolumeFace() = default;

    //! Constructor with intersection
    template <class Intersection, class GeometryHelper>
    FreeFlowStaggeredSubControlVolumeFace(const Intersection& is,
                               const typename Intersection::Geometry& isGeometry,
                               GridIndexType scvfIndex,
                               const std::vector<GridIndexType>& scvIndices,
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
      FreeFlowStaggeredSubControlVolumeFace(const GlobalPosition& dofPosition,
                                    const std::vector<GridIndexType>& scvIndices)
      {
          isGhostFace_ = true;
          center_ = dofPosition;
          scvIndices_ = scvIndices;
          scvfIndex_ = -1;
          dofIdx_ = -1;
      }

    //! The center of the sub control volume face
    const GlobalPosition& center() const
    {
        return center_;
    }

    //! The center of the sub control volume face
    const GlobalPosition& dofPosition() const
    {
        return center_;
    }

    //! The integration point for flux evaluations in global coordinates
    const GlobalPosition& ipGlobal() const
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

    const GlobalPosition& unitOuterNormal() const
    {
        return unitOuterNormal_;
    }

    //! index of the inside sub control volume for spatial param evaluation
    GridIndexType insideScvIdx() const
    {
        return scvIndices_[0];
    }

    //! index of the outside sub control volume for spatial param evaluation
    // This results in undefined behaviour if boundary is true
    GridIndexType outsideScvIdx() const
    {
        return scvIndices_[1];
    }

    //! The global index of this sub control volume face
    GridIndexType index() const
    {
        return scvfIndex_;
    }

    const GlobalPosition& corner(unsigned int localIdx) const
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
    GridIndexType dofIndex() const
    {
        return dofIdx_;
    }

    //! The global index of the dof living on the opposing face
    GridIndexType dofIndexOpposingFace() const
    {
        return oppositeIdx_;
    }

    //! The local index of this sub control volume face
    GridIndexType localFaceIdx() const
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

    //! Returns the sign of the unit outer normal's vector
    Scalar directionSign() const
    {
        return outerNormalScalar_;
    }


    const PairData<Scalar, GlobalPosition>& pairData(const int idx) const
    {
        return pairData_[idx];
    }

    const auto& pairData() const
    {
        return pairData_;
    }

private:
    Dune::GeometryType geomType_;
    std::vector<GlobalPosition> corners_;
    Scalar area_;
    GlobalPosition center_;
    GlobalPosition unitOuterNormal_;
    GridIndexType scvfIndex_;
    std::vector<GridIndexType> scvIndices_;
    bool boundary_;

    int dofIdx_;
    int oppositeIdx_;
    Scalar selfToOppositeDistance_;
    std::array<PairData<Scalar, GlobalPosition>, numPairs> pairData_;
    int localFaceIdx_;
    int dirIdx_;
    bool normalInPosCoordDir_;
    Scalar outerNormalScalar_;
    bool isGhostFace_;
};



} // end namespace

#endif
