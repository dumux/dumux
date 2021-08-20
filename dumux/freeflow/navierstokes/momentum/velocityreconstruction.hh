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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::StaggeredVelocityReconstruction
 */
#ifndef DUMUX_NAVIERSTOKES_STAGGERED_VELOCITYRECONSTRUCTION_HH
#define DUMUX_NAVIERSTOKES_STAGGERED_VELOCITYRECONSTRUCTION_HH

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
        static_assert(FVElementGeometry::GridGeometry::discMethod == DiscretizationMethod::cctpfa);
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
