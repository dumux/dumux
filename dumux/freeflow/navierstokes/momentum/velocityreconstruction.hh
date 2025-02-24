// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::StaggeredVelocityReconstruction
 */
#ifndef DUMUX_NAVIERSTOKES_STAGGERED_VELOCITYRECONSTRUCTION_HH
#define DUMUX_NAVIERSTOKES_STAGGERED_VELOCITYRECONSTRUCTION_HH

#include <numeric>
#include <algorithm>
#include <dune/common/reservedvector.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Helper class for reconstructing the velocity.
 */
struct StaggeredVelocityReconstruction
{
    //! Return the velocity vector at the center of the primal grid.
    template<class VelocityHelper, class FVElementGeometry>
    static auto cellCenterVelocity(const VelocityHelper& getFaceVelocity,
                                   const FVElementGeometry& fvGeometry)
    {
        static_assert(FVElementGeometry::GridGeometry::discMethod == DiscretizationMethods::cctpfa);
        using VelocityVector = typename FVElementGeometry::GridGeometry::GlobalCoordinate;
        VelocityVector result(0.0);

        const auto directionIndex = [&](const auto& vector)
        {
            return std::find_if(vector.begin(), vector.end(), [](const auto& x) { return std::abs(x) > 1e-8; } ) - vector.begin();
        };

        for (const auto& scvf : scvfs(fvGeometry))
        {
            const auto dirIdx = directionIndex(scvf.unitOuterNormal());
            result[dirIdx] += 0.5*getFaceVelocity(fvGeometry, scvf)[dirIdx];
        }
        return result;
    }

    //! Return the velocity vector at dof position given an scv
    template<class SubControlVolume, class FVElementGeometry,  class VelocityHelper>
    static auto faceVelocity(const SubControlVolume& scv,
                             const FVElementGeometry& fvGeometry,
                             const VelocityHelper& getNormalVelocityDofValue)
    {
        int axis = scv.dofAxis();
        const int dim = FVElementGeometry::GridGeometry::GridView::dimension;
        using Scalar = typename SubControlVolume::Traits::Scalar;
        using GlobalPosition = typename FVElementGeometry::GridGeometry::GlobalCoordinate;
        using VelocityVector = typename FVElementGeometry::GridGeometry::GlobalCoordinate;
        VelocityVector faceVelocityVector(0.0);

        // per dimension, we have at max two velocities from which we'll compute an average
        std::array<Dune::ReservedVector<Scalar, 2>, dim> normalVelocities;
        if (scv.boundary() && !fvGeometry.gridGeometry().dofOnPeriodicBoundary(scv.dofIndex()))
        {
            // iterate through the inner lateral velocities,
            for (const auto& scvf : scvfs(fvGeometry, scv))
            {
                if (scvf.isFrontal())
                    continue;

                // at a lateral velocity, find the inner and outer normal velocities
                const auto& orthogonalScvf = fvGeometry.lateralOrthogonalScvf(scvf);
                const auto& lateralScv = fvGeometry.scv(orthogonalScvf.insideScvIdx());
                auto lateralAxis = lateralScv.dofAxis();
                normalVelocities[lateralAxis].push_back( getNormalVelocityDofValue(lateralScv) ) ;
            }
        }
        else
        {
            // Find the location of interpolation, if periodic, from both sides
            const GlobalPosition selfPosition = scv.dofPosition();
            GlobalPosition outsideSelfPosition = selfPosition;
            if (fvGeometry.gridGeometry().dofOnPeriodicBoundary(scv.dofIndex()))
                outsideSelfPosition = fvGeometry.outsidePeriodicScv(scv).dofPosition();

            // iterate through the lateral velocities to get to the normal locations
            for (const auto& scvf : scvfs(fvGeometry, scv))
            {
                if (scvf.isFrontal())
                    continue;

                // at a lateral velocity, find the inner and outer normal velocities
                const auto& orthogonalScvf = fvGeometry.lateralOrthogonalScvf(scvf);
                const auto& insideLateralScv = fvGeometry.scv(orthogonalScvf.insideScvIdx());
                const auto& outsideLateralScv = fvGeometry.scv(orthogonalScvf.outsideScvIdx());
                const auto& lateralAxis = insideLateralScv.dofAxis();

                // Find the inside normal velocities
                const auto& insideNormalVelocity = getNormalVelocityDofValue(insideLateralScv);
                const auto& insideNormalPosition = insideLateralScv.dofPosition()[axis];

                // Find the outside normal velocities
                const auto& outsideNormalVelocity = getNormalVelocityDofValue(outsideLateralScv);
                const auto& outsideNormalPosition = outsideLateralScv.dofPosition()[axis];

                // Linear interpolation at the face plane and add to normal velocity collection
                const auto& innerDistance = std::abs(insideNormalPosition - selfPosition[axis]);
                const auto& totalDistance = std::abs(outsideNormalPosition - outsideSelfPosition[axis]) + innerDistance;
                const auto& velDiff = outsideNormalVelocity - insideNormalVelocity;

                normalVelocities[lateralAxis].push_back(insideNormalVelocity + (velDiff * innerDistance / totalDistance));
            }
        }

        for (int i = 0; i < faceVelocityVector.size(); i++)
        {
            if (i == axis)
                faceVelocityVector[i] = getNormalVelocityDofValue(scv);
            else
                faceVelocityVector[i] = std::accumulate(normalVelocities[i].begin(), normalVelocities[i].end(), 0.0) / normalVelocities[i].size();
        }

        return faceVelocityVector;
    }

};

} // end namespace Dumux

#endif // DUMUX_NAVIERSTOKES_STAGGERED_VELOCITYRECONSTRUCTION_HH
