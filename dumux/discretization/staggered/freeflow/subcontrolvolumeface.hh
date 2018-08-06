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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::FreeFlowStaggeredSubControlVolumeFace
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FREE_FLOW_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_STAGGERED_FREE_FLOW_SUBCONTROLVOLUMEFACE_HH

#include <utility>
#include <dune/common/fvector.hh>

#include <dumux/discretization/subcontrolvolumefacebase.hh>
#include <dumux/discretization/staggered/subcontrolvolumeface.hh>
#include <dumux/discretization/staggered/freeflow/staggeredgeometryhelper.hh>

#include <typeinfo>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Class for a sub control volume face in the staggered method, i.e a part of the boundary
 *        of a sub control volume we compute fluxes on. This is a specialization for free flow models.
 */
template<class GV,
         class T = StaggeredDefaultScvfGeometryTraits<GV> >
class FreeFlowStaggeredSubControlVolumeFace
: public SubControlVolumeFaceBase<FreeFlowStaggeredSubControlVolumeFace<GV, T>, T>
{
    using ThisType = FreeFlowStaggeredSubControlVolumeFace<GV, T>;
    using ParentType = SubControlVolumeFaceBase<ThisType, T>;
    using Geometry = typename T::Geometry;
    using GridIndexType = typename T::GridIndexType;

    using Scalar = typename T::Scalar;
    static const int dim = Geometry::mydimension;
    static const int dimworld = Geometry::coorddimension;

    static constexpr int numPairs = 2 * (dimworld - 1);

public:
    using GlobalPosition = typename T::GlobalPosition;

    //! State the traits public and thus export all types
    using Traits = T;

    // The default constructor
    FreeFlowStaggeredSubControlVolumeFace() = default;

    //! Constructor with intersection
    template <class Intersection, class GeometryHelper>
    FreeFlowStaggeredSubControlVolumeFace(const Intersection& is,
                                          const typename Intersection::Geometry& isGeometry,
                                          GridIndexType scvfIndex,
                                          const std::vector<GridIndexType>& scvIndices,
                                          const GeometryHelper& geometryHelper)
    : ParentType(),
      geomType_(isGeometry.type()),
      area_(isGeometry.volume()),
      center_(isGeometry.center()),
      unitOuterNormal_(is.centerUnitOuterNormal()),
      scvfIndex_(scvfIndex),
      scvIndices_(scvIndices),
      boundary_(is.boundary()),

      axisData_(geometryHelper.axisData()),
      pairData_(geometryHelper.pairData()),
      localFaceIdx_(geometryHelper.localFaceIndex()),
      dirIdx_(geometryHelper.directionIndex()),
      outerNormalSign_(sign(unitOuterNormal_[directionIndex()])),
      isGhostFace_(false)
      {
          corners_.resize(isGeometry.corners());
          for (int i = 0; i < isGeometry.corners(); ++i)
              corners_[i] = isGeometry.corner(i);
      }

    //! Constructor for a ghost face outside of the domain. Only needed to retrieve the center and scvIndices
    FreeFlowStaggeredSubControlVolumeFace(const GlobalPosition& dofPosition,
                                          const std::vector<GridIndexType>& scvIndices,
                                          const unsigned int dirIdx,
                                          const int dofIdx,
                                          const int scvfIndex)
    : center_(dofPosition),
      scvfIndex_(scvfIndex),
      scvIndices_(scvIndices),
      dofIdx_(dofIdx),
      selfToOppositeDistance_(0.0),
      dirIdx_(dirIdx),
      isGhostFace_(true)
    {}

    //! The center of the sub control volume face
    const GlobalPosition& center() const
    {
        return center_;
    }

    //! The position of the dof living on the face
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

    //! Returns bolean if the sub control volume face is on the boundary
    bool boundary() const
    {
        return boundary_;
    }

    //! The unit outer normal vector
    const GlobalPosition& unitOuterNormal() const
    {
        return unitOuterNormal_;
    }

    //! Index of the inside sub control volume for spatial param evaluation
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

    //! The positions of the corners
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

    //! The local index of this sub control volume face
    GridIndexType localFaceIdx() const
    {
        return localFaceIdx_;
    }

    //! Returns the dirction index of the facet (0 = x, 1 = y, 2 = z)
    unsigned int directionIndex() const
    {
        return dirIdx_;
    }

    //! Returns whether the unitNormal of the face points in positive coordinate direction
    bool normalInPosCoordDir() const
    {
        return directionSign() > 0;
    }

    //! Returns the sign of the unit outer normal's vector
    int directionSign() const
    {
        return outerNormalSign_;
    }

    //! Returns the data for one sub face
    const PairData<Scalar, GlobalPosition>& pairData(const int idx) const
    {
        return pairData_[idx];
    }

    //! Return an array of all pair data
    const auto& pairData() const
    {
        return pairData_;
    }

    //! Return an array of all pair data
    const auto& axisData() const
    {
        return axisData_;
    }

    bool isGhostFace() const
    {
        return isGhostFace_;
    }
    bool hasFirstParallelNeighbor(const int localSubFaceIdx) const
    {
        return !(pairData(localSubFaceIdx).parallelDofs[0] < 0);
    }

    bool hasSecondParallelNeighbor(const int localSubFaceIdx) const
    {
        return !(pairData(localSubFaceIdx).parallelDofs[1] < 0);
    }

    bool hasOuterNormal(const int localSubFaceIdx) const
    {
        return !(pairData_[localSubFaceIdx].normalPair.second < 0);
    }

    bool hasBackwardNeighbor() const
    {
      return  (axisData().inAxisBackwardDofs[0] >= 0);
    }

    bool hasForwardNeighbor() const
    {
        return (axisData().inAxisForwardDofs[0] >= 0);
    }

    bool canSecondOrder() const
    {
        if( hasBackwardNeighbor() && hasForwardNeighbor() && hasFirstParallelNeighbor() && hasSecondParallelNeighbor() )
            return 1;
        else
            return 0;
    }

private:

    bool hasOuterNormal() const
    {
        std::vector<Scalar> outerNormalDofIdxs;
        for (int i = 0; i < pairData().size() ; i++)
        {
            outerNormalDofIdxs.push_back(pairData_[i].normalPair.second);
        }
        if (std::any_of(outerNormalDofIdxs.begin(), outerNormalDofIdxs.end(), [](int j){return j < 0;}))
            return 0;
        else
            return 1;
    }

    bool hasFirstParallelNeighbor() const
    {
        std::vector<Scalar> firstParallelDofIdxs;
        for (int i = 0; i < pairData().size() ; i++)
        {
            firstParallelDofIdxs.push_back(pairData_[i].firstParallelFaceDofIdx);
        }
        if (std::any_of(firstParallelDofIdxs.begin(), firstParallelDofIdxs.end(), [](int j){return j < 0;}))
            return 0;
        else
            return 1;
    }

    bool hasSecondParallelNeighbor() const
    {
        std::vector<Scalar> secondParallelDofIdxs;
        for (int i = 0; i < pairData().size() ; i++)
        {
            secondParallelDofIdxs.push_back(pairData_[i].secondParallelFaceDofIdx);
        }
        if (std::any_of(secondParallelDofIdxs.begin(), secondParallelDofIdxs.end(), [](int j){return j < 0;}))
          return 0;
        else
          return 1;
    }

    Dune::GeometryType geomType_;
    std::vector<GlobalPosition> corners_;
    Scalar area_;
    GlobalPosition center_;
    GlobalPosition unitOuterNormal_;
    GridIndexType scvfIndex_;
    std::vector<GridIndexType> scvIndices_;
    bool boundary_;

    int dofIdx_;
    Scalar selfToOppositeDistance_;
    AxisData<Scalar> axisData_;
    std::array<PairData<Scalar, GlobalPosition>, numPairs> pairData_;

    int localFaceIdx_;
    unsigned int dirIdx_;
    int outerNormalSign_;
    bool isGhostFace_;
};

} // end namespace Dumux

#endif
