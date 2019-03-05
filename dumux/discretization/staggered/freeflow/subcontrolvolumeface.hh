// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/subcontrolvolumefacebase.hh>
#include <dumux/discretization/staggered/subcontrolvolumeface.hh>
#include <dumux/discretization/staggered/freeflow/staggeredgeometryhelper.hh>

#include <typeinfo>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Default traits class to be used for the sub-control volume faces
 *        for the free-flow staggered finite volume scheme
 * \tparam GridView the type of the grid view
 * \tparam upwindSchemeOrder the order of the upwind scheme
 */
template<class GridView, int upwindSchemeOrder>
struct FreeFlowStaggeredDefaultScvfGeometryTraits
{
    using Geometry = typename GridView::template Codim<1>::Geometry;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using Scalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;
    using PairData = typename FreeFlowStaggeredGeometryHelper<GridView, upwindSchemeOrder>::PairData;
    using AxisData = typename FreeFlowStaggeredGeometryHelper<GridView, upwindSchemeOrder>::AxisData;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Class for a sub control volume face in the staggered method, i.e a part of the boundary
 *        of a sub control volume we compute fluxes on. This is a specialization for free flow models.
 */
template<class GV,
         int upwindSchemeOrder,
         class T = FreeFlowStaggeredDefaultScvfGeometryTraits<GV, upwindSchemeOrder>>
class FreeFlowStaggeredSubControlVolumeFace
: public SubControlVolumeFaceBase<FreeFlowStaggeredSubControlVolumeFace<GV, upwindSchemeOrder, T>, T>
{
    using ThisType = FreeFlowStaggeredSubControlVolumeFace<GV, upwindSchemeOrder, T>;
    using ParentType = SubControlVolumeFaceBase<ThisType, T>;
    using Geometry = typename T::Geometry;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using LocalIndexType = typename IndexTraits<GV>::LocalIndex;

    using Scalar = typename T::Scalar;
    static const int dim = Geometry::mydimension;
    static const int dimworld = Geometry::coorddimension;

    static constexpr int numPairs = 2 * (dimworld - 1);

    static constexpr bool useHigherOrder = upwindSchemeOrder > 1;

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
      {    axisData_.selfDof = dofIdx; }

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
    LocalIndexType localFaceIdx() const
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
    const PairData<Scalar, GlobalPosition, upwindSchemeOrder>& pairData(const int idx) const
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

    //! Returns @c true if the face is a ghost face
    bool isGhostFace() const
    {
        return isGhostFace_;
    }

   /*!
    * \brief Check if the face has a parallel neighbor
    *
    * \param localSubFaceIdx The local index of the subface
    * \param parallelDegreeIdx The index describing how many faces away from the self face
    */
    bool hasParallelNeighbor(const int localSubFaceIdx, const int parallelDegreeIdx) const
    {
        return !(pairData(localSubFaceIdx).parallelDofs[parallelDegreeIdx] < 0);
    }

   /*!
    * \brief Check if the face has an outer normal neighbor
    *
    * \param localSubFaceIdx The local index of the subface
    */
    bool hasOuterNormal(const int localSubFaceIdx) const
    {
        return !(pairData_[localSubFaceIdx].normalPair.second < 0);
    }

   /*!
    * \brief Check if the face has a backward neighbor
    *
    * \param backwardIdx The index describing how many faces backward this dof is from the opposite face
    */
    template<bool enable = useHigherOrder, std::enable_if_t<enable, int> = 0>
    bool hasBackwardNeighbor(const int backwardIdx) const
    {
        return !(axisData().inAxisBackwardDofs[backwardIdx] < 0);
    }

   /*!
    * \brief Check if the face has a forward neighbor
    *
    * \param forwardIdx The index describing how many faces forward this dof is of the self face
    */
    template<bool enable = useHigherOrder, std::enable_if_t<enable, int> = 0>
    bool hasForwardNeighbor(const int forwardIdx) const
    {
        return !(axisData().inAxisForwardDofs[forwardIdx] < 0);
    }

    //! Returns the dof of the face
    GridIndexType dofIndex() const
    {
        return axisData().selfDof;
    }

    //! Returns the dof of the opposing face
    GridIndexType dofIndexOpposingFace() const
    {
        return axisData().oppositeDof;
    }

    //! Returns the dof the first forward face
    GridIndexType dofIndexForwardFace() const
    {
        return axisData().inAxisForwardDofs[0];
    }

    //! Returns the dof of the first backward face
    GridIndexType dofIndexBackwardFace() const
    {
        return axisData().inAxisBackwardDofs[0];
    }

    //! Returns the distance between the face and the opposite one
    Scalar selfToOppositeDistance() const
    {
        return axisData().selfToOppositeDistance;
    }

    /*!
    * \brief Returns the distance between the parallel dofs
    *
    * \param localSubFaceIdx The local index of the subface
    * \param parallelDegreeIdx The index describing how many faces away from the self
    */
    Scalar cellCenteredParallelDistance(const int localSubFaceIdx, const int parallelDegreeIdx) const
    {
        return (pairData(localSubFaceIdx).parallelDistances[parallelDegreeIdx] +
                pairData(localSubFaceIdx).parallelDistances[parallelDegreeIdx+1]) * 0.5;
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
    Scalar selfToOppositeDistance_;
    AxisData<Scalar, upwindSchemeOrder> axisData_;
    std::array<PairData<Scalar, GlobalPosition, upwindSchemeOrder>, numPairs> pairData_;

    int localFaceIdx_;
    unsigned int dirIdx_;
    int outerNormalSign_;
    bool isGhostFace_;
};

} // end namespace Dumux

#endif // DUMUX_DISCRETIZATION_STAGGERED_FREE_FLOW_SUBCONTROLVOLUMEFACE_HH
