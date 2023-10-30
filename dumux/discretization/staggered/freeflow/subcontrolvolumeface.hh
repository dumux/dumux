// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::FreeFlowStaggeredSubControlVolumeFace
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FREE_FLOW_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_STAGGERED_FREE_FLOW_SUBCONTROLVOLUMEFACE_HH

#include <array>
#include <utility>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

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
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using Scalar = typename GridView::ctype;
    using GeometryHelper = FreeFlowStaggeredGeometryHelper<GridView, upwindSchemeOrder>;
    using PairData = typename GeometryHelper::PairData;
    using AxisData = typename GeometryHelper::AxisData;

    using Grid = typename GridView::Grid;
    static constexpr int dim = Grid::dimension;
    static constexpr int dimWorld = Grid::dimensionworld;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Geometry = Dune::AxisAlignedCubeGeometry<Scalar, dim-1, dimWorld>;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Helper function to turn a given cell scvface into a fake boundary face
 * \note This function is considered internal to staggered freeflow and may change or be deleted at any time without deprecation warning
 */
template<class SubControlVolumeFace>
SubControlVolumeFace makeStaggeredBoundaryFace(const SubControlVolumeFace& scvf,
                                               const typename SubControlVolumeFace::GlobalPosition& newCenter)
{
    auto bf = scvf;
    bf.setCenter(newCenter);
    bf.setBoundary(true);
    bf.setIsGhostFace(true);
    return bf;
}

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

    using PairData = typename T::PairData;
    using AxisData = typename T::AxisData;

    using Scalar = typename T::Scalar;
    static const int dim = GV::dimension;

    static constexpr int numPairs = 2 * (dim - 1);

    static constexpr bool useHigherOrder = upwindSchemeOrder > 1;

public:
    using GlobalPosition = typename T::GlobalPosition;

    static constexpr int numCornersPerFace = 2 * (dim - 1);

    //! State the traits public and thus export all types
    using Traits = T;

    // The default constructor
    FreeFlowStaggeredSubControlVolumeFace() = default;

    //! Constructor with intersection
    template <class Intersection>
    FreeFlowStaggeredSubControlVolumeFace(const Intersection& is,
                                          const typename Intersection::Geometry& isGeometry,
                                          GridIndexType scvfIndex,
                                          const std::vector<GridIndexType>& scvIndices,
                                          const typename T::GeometryHelper& geometryHelper)
    : ParentType(),
      area_(isGeometry.volume()),
      center_(isGeometry.center()),
      unitOuterNormal_(is.centerUnitOuterNormal()),
      scvfIndex_(scvfIndex),
      scvIndices_(scvIndices),
      boundary_(is.boundary()),
      axisData_(geometryHelper.axisData()),
      pairData_(std::move(geometryHelper.pairData())),
      localFaceIdx_(geometryHelper.localFaceIndex()),
      dirIdx_(geometryHelper.directionIndex()),
      outerNormalSign_(sign(unitOuterNormal_[directionIndex()])),
      isGhostFace_(false)
    {
        dimensions[0] = (isGeometry.corner(1) - isGeometry.corner(0)).two_norm();
        if constexpr (dim == 3)
            dimensions[1] = (isGeometry.corner(2) - isGeometry.corner(0)).two_norm();
    }

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

    //! Returns boolean if the sub control volume face is on the boundary
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
    GridIndexType outsideScvIdx() const
    {
        return scvIndices_[1];
    }

    //! The global index of this sub control volume face
    GridIndexType index() const
    {
        return scvfIndex_;
    }

    //! The local index of this sub control volume face
    LocalIndexType localFaceIdx() const
    {
        return localFaceIdx_;
    }

    //! Returns the direction index of the facet (0 = x, 1 = y, 2 = z)
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
    const PairData& pairData(const int idx) const
    {
        return pairData_[idx];
    }

    //! Return an array of all pair data
    const std::array<PairData, numPairs>& pairData() const
    {
        return pairData_;
    }

    //! Return an array of all pair data
    const AxisData& axisData() const
    {
        return axisData_;
    }

    //! Returns true if the face is a ghost face
    bool isGhostFace() const
    {
        return isGhostFace_;
    }

    //! Returns the length of the face in a certain direction (adaptation of area() for 3d)
    Scalar faceLength(const int localSubFaceIdx) const
    {
        if (dim == 3 && localSubFaceIdx > 1)
            return dimensions[1];
        else
            return dimensions[0];
    }

   /*!
    * \brief Check if the face has a parallel neighbor
    *
    * \param localSubFaceIdx The local index of the subface
    * \param parallelDegreeIdx The index describing how many faces away from the self face
    */
    bool hasParallelNeighbor(const int localSubFaceIdx, const int parallelDegreeIdx) const
    {
        return pairData(localSubFaceIdx).hasParallelNeighbor[parallelDegreeIdx];
    }

    /*!
    * \brief Check if the face has a half parallel neighbor
    *
    * \param localSubFaceIdx The local index of the subface
    *
    * ------------
    * |          |
    * |          |
    * |          |
    * -----------------------
    * | yyyyyyyy s          |
    * | yyyyyyyy s          |
    * | yyyyyyyy s          |
    * -----------------------
    * In this corner geometry, hasParallelNeighbor will return true for subcontrolvolumeface s belonging to the
    * element filled by 'y's, but hasParallelNeighbor will return false for the subcontrolvolumeface that has the
    * same dofIndex. We name this situation hasHalfParallelNeighbor.
    */
    bool hasHalfParallelNeighbor(const int localSubFaceIdx) const
    {
        return pairData(localSubFaceIdx).hasHalfParallelNeighbor;
    }

    /*!
    * \brief Check if the face has a corner parallel neighbor
    *
    * \param localSubFaceIdx The local index of the subface
    *
    * ------------
    * | yyyyyyyy s
    * | yyyyyyyy s
    * | yyyyyyyy s
    * -----------------------
    * |          |          |
    * |          |          |
    * |          |          |
    * -----------------------
    * In this corner geometry, hasParallelNeighbor will return true for subcontrolvolumeface s belonging to the
    * element filled by 'y's. However, as there also might be a boundary velocity value known at the corner, which
    * can be used instead of the standard parallel velocity in some cases, we want to identify this situation. We
    * name it cornerParallelNeighbor.
    */
    bool hasCornerParallelNeighbor(const int localSubFaceIdx) const
    {
        return pairData(localSubFaceIdx).hasCornerParallelNeighbor;
    }

   /*!
    * \brief Check if the face has an outer normal neighbor
    *
    * \param localSubFaceIdx The local index of the subface
    */
    bool hasOuterLateral(const int localSubFaceIdx) const
    {
        return pairData(localSubFaceIdx).hasOuterLateral;
    }

   /*!
    * \brief Check if the face has a backward neighbor
    *
    * \param backwardIdx The index describing how many faces backward this dof is from the opposite face
    */
    template<bool enable = useHigherOrder, std::enable_if_t<enable, int> = 0>
    bool hasBackwardNeighbor(const int backwardIdx) const
    {
        return axisData().hasBackwardNeighbor[backwardIdx];
    }

   /*!
    * \brief Check if the face has a forward neighbor
    *
    * \param forwardIdx The index describing how many faces forward this dof is of the self face
    */
    template<bool enable = useHigherOrder, std::enable_if_t<enable, int> = 0>
    bool hasForwardNeighbor(const int forwardIdx) const
    {
        return axisData().hasForwardNeighbor[forwardIdx];
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
    Scalar parallelDofsDistance(const int localSubFaceIdx, const int parallelDegreeIdx) const
    {
        if (parallelDegreeIdx == 0)
            return (faceLength(localSubFaceIdx) + pairData(localSubFaceIdx).parallelCellWidths[0]) * 0.5;
            // pairData(localSubFaceIdx).parallelCellWidths[0]) will return 0.0 if the subface perpendicular the scvf lies on a boundary
        else
        {
            assert((parallelDegreeIdx == 1) && "Only the width of the first two parallel cells (indices 0 and 1) is stored for each scvf.");
            return (pairData(localSubFaceIdx).parallelCellWidths[0] + pairData(localSubFaceIdx).parallelCellWidths[1]) * 0.5;
        }
    }

    //! set the center to a different position
    void setCenter(const GlobalPosition& center)
    { center_ = center; }

    //! set the boundary flag
    void setBoundary(bool boundaryFlag)
    { boundary_ = boundaryFlag; }

    //! set the ghost face flag
    void setIsGhostFace(bool isGhostFaceFlag)
    { isGhostFace_ = isGhostFaceFlag; }

private:
    std::array<Scalar, dim-1> dimensions;
    Scalar area_;
    GlobalPosition center_;
    GlobalPosition unitOuterNormal_;
    GridIndexType scvfIndex_;
    std::vector<GridIndexType> scvIndices_;
    bool boundary_;

    Scalar selfToOppositeDistance_;
    AxisData axisData_;
    std::array<PairData, numPairs> pairData_;

    int localFaceIdx_;
    unsigned int dirIdx_;
    int outerNormalSign_;
    bool isGhostFace_;
};

} // end namespace Dumux

#endif
