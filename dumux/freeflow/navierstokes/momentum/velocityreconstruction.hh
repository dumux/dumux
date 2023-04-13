// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::StaggeredVelocityReconstruction
 */
#ifndef DUMUX_NAVIERSTOKES_STAGGERED_VELOCITYRECONSTRUCTION_HH
#define DUMUX_NAVIERSTOKES_STAGGERED_VELOCITYRECONSTRUCTION_HH
#include <algorithm>
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
};

} // end namespace Dumux

#endif // DUMUX_NAVIERSTOKES_STAGGERED_VELOCITYRECONSTRUCTION_HH
