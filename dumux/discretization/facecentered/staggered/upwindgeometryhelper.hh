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
 * \ingroup FaceCenteredStaggeredDiscretization
 * \copydoc Dumux::FaceCenteredStaggeredUpwindGeometryHelper
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_UPWIND_GEOMETRY_HELPER_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_UPWIND_GEOMETRY_HELPER_HH

#include <bitset>
#include <type_traits>
#include <utility>

#include <dune/common/reservedvector.hh>

#include <dumux/common/indextraits.hh>
#include <dune/common/iteratorrange.hh>
#include <dumux/discretization/scvandscvfiterators.hh>
#include <dumux/discretization/facecentered/staggered/normalaxis.hh>

namespace Dumux {


/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Geometry helper for (higher order) upwinding
 */
template<class FVElementGeometry>
class FaceCenteredStaggeredUpwindGeometryHelper
{
    using GG = typename FVElementGeometry::GridGeometry;
    using GridView = typename GG::GridView;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;

public:
    //! export type of subcontrol volume face
    using SubControlVolume = typename GG::SubControlVolume;
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    using Scalar = typename SubControlVolume::Traits::Scalar;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GridGeometry = GG;
    using UpwindScheme = typename GridGeometry::UpwindScheme;

    FaceCenteredStaggeredUpwindGeometryHelper(const FVElementGeometry& fvGeometry)
    : fvGeometry_(fvGeometry)
    {}

    const UpwindScheme& upwindMethods() const
    { return gridGeometry().upwindMethods(); }


    //! Return the scv in the neighbor element with the same local index
    const SubControlVolume& nextCorrespondingScv(const SubControlVolume& startScv, const LocalIndexType& targetLocalScvIndex ) const
    {
        assert(!(startScv.boundary()));
        auto nextFVGeometry = localView(gridGeometry());
        const auto& nextElement = gridGeometry().element(startScv.neighborElementIdx());
        nextFVGeometry.bind(nextElement);
        for (auto&& neighborElementScv : scvs(nextFVGeometry))
        {
            if (neighborElementScv.localDofIndex() == targetLocalScvIndex)
                return neighborElementScv;
        }
        DUNE_THROW(Dune::InvalidStateException, "No scv found");
    }


    //! Return a the lateral scvf opposite to the given sub control volume face
    const SubControlVolumeFace& oppositeLateralScvf(const SubControlVolumeFace& lateralScvf) const
    {
        assert(lateralScvf.isLateral());
        const auto& selfScv = scv(lateralScvf.insideScvIdx());
        GridIndexType oppositeLateralScvfIdx = -1;

        for (auto&& scvf : scvfs(*this, selfScv))
        {
            if (scvf.isLateral() // must be a lateral face
                && scvf.localIndex() != lateralScvf.localIndex() // should not be the same face
                && scvf.normalAxis() == lateralScvf.normalAxis() /* should have the same normal direction index*/)
            {
                oppositeLateralScvfIdx = scvf.index();
            }
        }
        return scvf(oppositeLateralScvfIdx);
    }

    //! Returns the lateral length of the scv in a certain direction
    Scalar scvLateralLength(const SubControlVolume& scv, const LocalIndexType& normalAxis) const
    {
        auto fvGeometry = localView(gridGeometry());
        const auto& element = fvGeometry.gridGeometry().element(scv.elementIndex());
        fvGeometry.bind(element);

        std::array<SubControlVolumeFace,2> opposingLateralFaces;
        LocalIndexType i = 0;
        for (auto&& scvf : scvfs(fvGeometry, scv))
        {
            if (scvf.isLateral() && scvf.normalAxis() == normalAxis ) // find a lateral face on the correct axis (there are 2)
            {
                opposingLateralFaces[i] = scvf;
                i++;
            }
        }
        return (opposingLateralFaces[0].center() - opposingLateralFaces[1].center()).two_norm();
    }

    //! Returns the frontal length of the scv (half of selfToOppositeDistance)
    Scalar scvFrontalLength(const SubControlVolume& scv) const
    {
        auto fvGeometry = localView(gridGeometry());
        const auto& element = fvGeometry.gridGeometry().element(scv.elementIndex());
        fvGeometry.bind(element);

        GlobalPosition frontalFaceLocation(0.0);
        for (auto&& scvf : scvfs(fvGeometry, scv))
        {
            if (scvf.isFrontal() && !scvf.boundary())
                frontalFaceLocation = scvf.center();
        }

        return (scv.dofPosition() - frontalFaceLocation).two_norm();
    }

    /////////////////////////
    /// forward neighbors ///
    /////////////////////////

    bool hasForwardNeighbor(const SubControlVolumeFace& frontalScvf) const
    {
        assert(frontalScvf.isFrontal());
        const auto& selfScv = scv(frontalScvf.insideScvIdx());
        return !selfScv.boundary();
    }

    GridIndexType forwardScvIdx(const SubControlVolumeFace& frontalScvf) const
    {
        assert(frontalScvf.isFrontal());
        assert(hasForwardNeighbor(frontalScvf));
        const auto& selfScv = scv(frontalScvf.insideScvIdx());
        const auto& nextScv = nextCorrespondingScv(selfScv, selfScv.indexInElement());
        return nextScv.index();
    }

    //////////////////////////
    /// backward neighbors ///
    //////////////////////////

    bool hasBackwardNeighbor(const SubControlVolumeFace& frontalScvf) const
    {
        assert(frontalScvf.isFrontal());
        const auto& oppositeScv = scv(frontalScvf.outsideScvIdx());
        return !oppositeScv.boundary();
    }

    GridIndexType backwardScvIdx(const SubControlVolumeFace& frontalScvf) const
    {
        assert(frontalScvf.isFrontal());
        assert(hasBackwardNeighbor(frontalScvf));
        const auto& oppositeScv = scv(frontalScvf.outsideScvIdx());
        const auto& nextScv = nextCorrespondingScv(oppositeScv, oppositeScv.indexInElement());
        return nextScv.index();
    }

    /////////////////////////
    /// Frontal Distances ///
    /////////////////////////

    Scalar selfToOppositeDistance(const SubControlVolumeFace& frontalScvf) const
    {
        assert(frontalScvf.isFrontal());
        const auto& selfScv = fvGeometry_.scv(frontalScvf.insideScvIdx());
        const auto& outsideScv = fvGeometry_.scv(frontalScvf.outsideScvIdx());
        return (selfScv.dofPosition() - outsideScv.dofPosition()).two_norm();
    }

    Scalar selfToForwardDistance(const SubControlVolumeFace& frontalScvf) const
    {
        assert(frontalScvf.isFrontal());

        const auto& selfScv = fvGeometry_.scv(frontalScvf.insideScvIdx());
        if (selfScv.boundary())
            return 0.0;

        const auto& forwardScv = nextCorrespondingScv(selfScv, selfScv.indexInElement());
        return (selfScv.dofPosition() - forwardScv.dofPosition()).two_norm();
    }

    Scalar oppositeToBackwardDistance(const SubControlVolumeFace& frontalScvf) const
    {
        assert(frontalScvf.isFrontal());
        const auto& outsideScv = fvGeometry_.scv(frontalScvf.outsideScvIdx());
        if (outsideScv.boundary())
            return 0.0;

        const auto& backwardScv = nextCorrespondingScv(outsideScv, outsideScv.indexInElement());
        return (outsideScv.dofPosition() - backwardScv.dofPosition()).two_norm();
    }

    //////////////////////////
    /// Parallel Neighbors ///
    //////////////////////////

    bool hasParallelNeighbor(const SubControlVolumeFace& lateralScvf) const
    {
        assert(lateralScvf.isLateral());
        return !lateralScvf.boundary();
    }

    bool hasSecondParallelNeighbor(const SubControlVolumeFace& lateralScvf) const
    {
        assert(lateralScvf.isLateral());
        assert(hasParallelNeighbor(lateralScvf));

        const auto& orthogonalScvf = lateralOrthogonalScvf(lateralScvf);
        const auto& orthagonalScv = scv(orthogonalScvf.insideScvIdx());
        const auto& nextOrthagonalScv = nextCorrespondingScv(orthagonalScv, orthagonalScv.indexInElement());
        return !nextOrthagonalScv.boundary();
    }

    GridIndexType parallelScvIdx(const SubControlVolumeFace& lateralScvf) const
    {
        assert(lateralScvf.isLateral());
        assert(hasParallelNeighbor(lateralScvf));
        return lateralScvf.outsideScvIdx();
    }

    GridIndexType secondParallelScvIdx(const SubControlVolumeFace& lateralScvf) const
    {
        assert(lateralScvf.isLateral());
        assert(hasSecondParallelNeighbor(lateralScvf));

        const auto& selfScv = scv(lateralScvf.insideScvIdx());
        const auto& orthogonalScvf = lateralOrthogonalScvf(lateralScvf);
        const auto& orthagonalScv = scv(orthogonalScvf.insideScvIdx());
        const auto& nextOrthagonalScv = nextCorrespondingScv(orthagonalScv, orthagonalScv.indexInElement());
        const auto& secondParallelScv = nextCorrespondingScv(nextOrthagonalScv, selfScv.indexInElement());
        return secondParallelScv.index();
    }

    SubControlVolumeFace outerParallelLateralScvf(const SubControlVolumeFace& lateralScvf) const
    {
        assert(lateralScvf.isLateral());
        assert(hasParallelNeighbor(lateralScvf));

        const auto& parallelScv = scv(lateralScvf.outsideScvIdx());

        auto fvGeometry = localView(gridGeometry());
        const auto& element = fvGeometry.gridGeometry().element(parallelScv.elementIndex());
        fvGeometry.bind(element);
        GridIndexType index = -1;
        for (auto&& scvf : scvfs(fvGeometry, parallelScv))
        {
            if (scvf.localIndex() == lateralScvf.localIndex())
                index = scvf.index();
        }

        return scvf(index);
    }

    //////////////////////////
    /// Parallel Distances ///
    //////////////////////////

    Scalar selfToParallelDistance(const SubControlVolumeFace& lateralScvf) const
    {
        assert(lateralScvf.isLateral());

        // if there is no parallel neighbor, the distance from the dof to the boundary is used
        if (!hasParallelNeighbor(lateralScvf))
            return insideScvLateralLength(lateralScvf) / 2.0;

        const auto& selfScv = fvGeometry_.scv(lateralScvf.insideScvIdx());
        const auto& parallelScv = fvGeometry_.scv(lateralScvf.outsideScvIdx());
        return (parallelScv.dofPosition() - selfScv.dofPosition()).two_norm();
    }

    Scalar paralleltoSecondParallelDistance(const SubControlVolumeFace& lateralScvf) const
    {
        assert(lateralScvf.isLateral());
        assert(hasParallelNeighbor(lateralScvf));

        // if there is no second parallel neighbor, the distance from the parallel dof to the boundary is used
        if (!hasSecondParallelNeighbor(lateralScvf))
            return outsideScvLateralLength(lateralScvf) / 2.0;

        const auto& parallelScv = scv(parallelScvIdx(lateralScvf));
        const auto& secondParallelScv = scv(secondParallelScvIdx(lateralScvf));
        return (secondParallelScv.dofPosition() - parallelScv.dofPosition()).two_norm();
    }

    Scalar insideScvLateralLength(const SubControlVolumeFace& lateralScvf) const
    {
        assert(lateralScvf.isLateral());
        const auto& selfScv = fvGeometry_.scv(lateralScvf.insideScvIdx());
        return scvLateralLength(selfScv, lateralScvf.normalAxis());
    }

    Scalar outsideScvLateralLength(const SubControlVolumeFace& lateralScvf) const
    {
        assert(lateralScvf.isLateral());
        assert(hasParallelNeighbor(lateralScvf));

        const auto& outsideScv = scv(lateralScvf.outsideScvIdx());
        return scvLateralLength(outsideScv, lateralScvf.normalAxis());
    }

    /////////////////////////
    /// Half/Corner Bools ///
    /////////////////////////

    /*     "Half Parallel Neighbor"                                  "Corner Parallel Neighbor"
     *
     *   -----------                                              --------------------
     *   |    pppp b                                              |    pppp p        |
     *   |    pppp b                s: selfScv,                   |    pppp p        |        s: selfScv,
     *   |    pppp b                L: lateralScvf                |    pppp p        |        L: lateralScvf
     *   -----LLLL * bbbbbbbb       p: parallelScv                -----LLLL * bbbbbbbb        p: parallelScv
     *   |    ssss s        |       b: domain boundary            |    ssss b                 b: domain boundary
     *   |    ssss s        |                                     |    ssss b
     *   |    ssss s        |                                     |    ssss b
     *   --------------------                                     -----------
     */

    bool hasHalfParallelNeighbor(const SubControlVolumeFace& lateralScvf) const
    {
        assert(lateralScvf.isLateral());
        const auto& selfScv = scv(lateralScvf.insideScvIdx());
        if (!lateralScvf.boundary() && !selfScv.boundary())
        {
            const auto& parallelScv = scv(parallelScvIdx(lateralScvf));
            if (parallelScv.boundary())
                return true;
        }

        return false;
    }

    bool hasCornerParallelNeighbor(const SubControlVolumeFace& lateralScvf) const
    {
        assert(lateralScvf.isLateral());
        const auto& selfScv = scv(lateralScvf.insideScvIdx());
        if (!lateralScvf.boundary() && selfScv.boundary())
        {
            const auto& parallelScv = scv(parallelScvIdx(lateralScvf));
            if (!parallelScv.boundary())
                return true;
        }

        return false;
    }

    bool lateralFaceContactsBoundary(const SubControlVolumeFace& lateralScvf) const
    {
        return ( lateralScvf.boundary() ||
                 hasHalfParallelNeighbor(lateralScvf) ||
                 hasCornerParallelNeighbor(lateralScvf) );
    }

    GlobalPosition cornerBoundaryPosition(const SubControlVolumeFace& lateralScvf) const
    {
        assert(lateralScvf.isLateral());
        assert(hasCornerParallelNeighbor(lateralScvf) || hasHalfParallelNeighbor(lateralScvf) || lateralScvf.boundary());

        const auto& selfScv = scv(lateralScvf.insideScvIdx());
        const Scalar halfFaceLength = scvFrontalLength(selfScv) / 2.0;
        const GlobalPosition shift = selfScv.innerUnitNormal() * halfFaceLength;
        return lateralScvf.center() + shift;
    }

    //! The grid geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    { return fvGeometry_.gridGeometry(); }


    //! The bound element
    const Element& element() const
    { return *fvGeometry_.element(); }

private:
    GridIndexType eIdx_;
    const FVElementGeometry& fvGeometry_;
};



} // end namespace Dumux

#endif
