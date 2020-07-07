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
 * \copydoc Dumux::StaggeredVelocityGradients
 */
#ifndef DUMUX_NAVIERSTOKES_STAGGERED_VELOCITYGRADIENTS_HH
#define DUMUX_NAVIERSTOKES_STAGGERED_VELOCITYGRADIENTS_HH

#include <optional>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Helper class for calculating the velocity gradients for the Navier-Stokes model using the staggered grid discretization.
 */
class StaggeredVelocityGradients
{

public:

    /*!
     * \brief Returns the in-axis velocity gradient.
     *
     * \verbatim
     *              ---------=======                 == and # staggered half-control-volume
     *              |       #      | current scv
     *              |       #      |                 # staggered face over which fluxes are calculated
     *   vel.Opp <~~|       O~~>   x~~~~> vel.Self
     *              |       #      |                 x dof position
     *              |       #      |
     *              --------========                 -- element
     *
     *                                               O position at which gradient is evaluated (integration point)
     * \endverbatim
     */
    template<class FVElementGeometry, class ElemVolVars>
    static auto velocityGradII(const FVElementGeometry fvGeometry,
                               const typename FVElementGeometry::SubControlVolumeFace& scvf,
                               const ElemVolVars& elemVolVars)
    {
        assert(scvf.isFrontal());
        // The velocities of the dof at interest and the one of the opposite scvf.
        const auto velocitySelf = elemVolVars[scvf.insideScvIdx()].velocity();
        const auto velocityOpposite = elemVolVars[scvf.outsideScvIdx()].velocity();
        const auto distance = (fvGeometry.scv(scvf.outsideScvIdx()).dofPosition() - fvGeometry.scv(scvf.insideScvIdx()).dofPosition()).two_norm();

        return (velocityOpposite - velocitySelf) / distance * scvf.directionSign();
    }

    /*!
     * \brief Returns the velocity gradient perpendicular to the orientation of our current scvf.
     *
     * \verbatim
     *              ----------------
     *              |              |outer                    || and # staggered half-control-volume (own element)
     *              |              |vel.        gradient
     *              |              |~~~~>       ------->     :: staggered half-control-volume (neighbor element)
     *              |              |             ------>
     *              | lateral scvf |              ----->     x dof position
     *              ---------######O:::::::::      ---->
     *              |      ||      |       ::       --->     -- elements
     *              |      ||      |       ::        -->
     *              |      || scv  x~~~~>  ::         ->     O position at which gradient is evaluated (integration point)
     *              |      ||      | inner ::
     *              |      ||      | vel.  ::
     *              ---------#######:::::::::
     *
     *
     * \endverbatim
     */
    template<class FVElementGeometry, class ElemVolVars>
    static auto velocityGradIJ(const FVElementGeometry fvGeometry,
                               const typename FVElementGeometry::SubControlVolumeFace& scvf,
                               const ElemVolVars& elemVolVars)
    {
        assert(scvf.isLateral());

        const auto innerVelocity = elemVolVars[scvf.insideScvIdx()].velocity();
        const auto outerVelocity = elemVolVars[scvf.outsideScvIdx()].velocity();

        const auto distance = scvf.boundary() ? (fvGeometry.scv(scvf.insideScvIdx()).dofPosition() - scvf.ipGlobal()).two_norm()
                                              : (fvGeometry.scv(scvf.insideScvIdx()).dofPosition() - fvGeometry.scv(scvf.outsideScvIdx()).dofPosition()).two_norm();

        return (outerVelocity - innerVelocity) / distance * scvf.directionSign();
    }

    /*!
     * \brief Returns the velocity gradient in line with our current scvf.
     *
     * \verbatim
     *                      ^       gradient
     *                      |  ^
     *                      |  |  ^
     *                      |  |  |  ^
     *                      |  |  |  |  ^
     *                      |  |  |  |  |  ^       || and # staggered half-control-volume (own element)
     *                      |  |  |  |  |  |
     *                                             :: staggered half-control-volume (neighbor element)
     *              ----------------
     *              |     inner    |      outer    x dof position (of own scv)
     *              |      vel.    |       vel.
     *              |       ^      |        ^      -- elements
     *              |       | lat. |        |
     *              |       | scvf |        |      O position at which gradient is evaluated (integration point)
     *              ---------######O:::::::::
     *              |      ||      |       ::
     *              |      ||      |       ::
     *              |      || scv  x       ::
     *              |      ||      |       ::
     *              |      ||      |       ::
     *              ---------#######:::::::::
     *
     *
     * \endverbatim
     */
    template<class FVElementGeometry, class ElemVolVars>
    static auto velocityGradJI(const FVElementGeometry fvGeometry,
                               const typename FVElementGeometry::SubControlVolumeFace& scvf,
                               const ElemVolVars& elemVolVars)
    {
        assert(scvf.isLateral());
        const auto& orthogonalScvf = fvGeometry.scvfWithCommonEntity(scvf);

        const auto innerVelocity = elemVolVars[orthogonalScvf.insideScvIdx()].velocity();
        const auto outerVelocity = elemVolVars[orthogonalScvf.outsideScvIdx()].velocity();

        const auto distance = orthogonalScvf.boundary() ? (fvGeometry.scv(orthogonalScvf.insideScvIdx()).dofPosition() - orthogonalScvf.ipGlobal()).two_norm()
                                                        : (fvGeometry.scv(orthogonalScvf.insideScvIdx()).dofPosition() - fvGeometry.scv(orthogonalScvf.outsideScvIdx()).dofPosition()).two_norm();

        return (outerVelocity - innerVelocity) / distance * orthogonalScvf.directionSign();
    }
};

} // end namespace Dumux

#endif // DUMUX_NAVIERSTOKES_STAGGERED_VELOCITYGRADIENTS_HH
