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
template<class Scalar, class GridGeometry, class BoundaryTypes, class Indices>
class StaggeredVelocityGradients
{
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:

    /*!
     * \brief Returns the in-axis velocity gradient.
     *
     * \verbatim
     *              ---------=======                 == and # staggered half-control-volume
     *              |       #      | current scvf
     *              |       #      |                 # staggered face over which fluxes are calculated
     *   vel.Opp <~~|       O~~>   x~~~~> vel.Self
     *              |       #      |                 x dof position
     *        scvf  |       #      |
     *              --------========                 -- element
     *
     *                                               O position at which gradient is evaluated
     * \endverbatim
     */
    template<class StaggeredFVElementGeometry, class FaceVariables>
    static Scalar velocityGradII(const StaggeredFVElementGeometry fvGeometry,
                                 const typename StaggeredFVElementGeometry::StaggeredSubControlVolumeFace& staggeredScvf,
                                 const FaceVariables& faceVars)
    {
        // The velocities of the dof at interest and the one of the opposite scvf.
        const Scalar velocitySelf = faceVars.velocitySelf();
        const Scalar velocityOpposite = faceVars.velocityOpposite();

        return ((velocityOpposite - velocitySelf) / fvGeometry.normalDistanceForGradient(staggeredScvf)) * staggeredScvf.directionSign();
    }

    /*!
     * \brief Returns the velocity gradient perpendicular to the orientation of our current scvf.
     *
     * \verbatim
     *              ----------------
     *              |              |vel.
     *              |              |Parallel
     *              |              |~~~~>       ------->
     *              |              |             ------>     * gradient
     *              |              |              ----->
     *       scvf   ---------######O:::::::::      ---->     || and # staggered half-control-volume (own element)
     *              |      ||      | curr. ::       --->
     *              |      ||      | scvf  ::        -->     :: staggered half-control-volume (neighbor element)
     *              |      ||      x~~~~>  ::         ->
     *              |      ||      | vel.  ::                # lateral staggered faces over which fluxes are calculated
     *        scvf  |      ||      | Self  ::
     *              ---------#######:::::::::                x dof position
     *                 scvf
     *                                                       -- elements
     *
     *                                                       O position at which gradient is evaluated
     * \endverbatim
     */
    template<class Problem, class StaggeredFVElementGeometry, class FaceVariables>
    static Scalar velocityGradIJ(const Problem& problem,
                                 const Element& element,
                                 const StaggeredFVElementGeometry staggeredFVGeometry,
                                 const typename StaggeredFVElementGeometry::StaggeredSubControlVolumeFace& staggeredLateralScvf,
                                 const FaceVariables& faceVars,
                                 const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
                                 const std::optional<BoundaryTypes>& lateralFaceBoundaryTypes)
    {
        // For the velocityGrad_ij derivative, get the velocities at the current (own) scvf
        // and at the parallel one at the neighboring scvf.
        const Scalar innerParallelVelocity = faceVars.velocitySelf();

        const auto& staggeredScv = staggeredFVGeometry.scv(staggeredLateralScvf.insideScvIdx());
        const auto directionIndex = staggeredScv.directionIndex();
        const auto localLateralFaceIdx = staggeredFVGeometry.localLateralFaceIndex(staggeredLateralScvf); // goes from 0 to dim-1

        const auto outerParallelVelocity = [&]()
        {
            if (!staggeredLateralScvf.boundary())
                return faceVars.velocityParallel(localLateralFaceIdx, 0);
            else if (lateralFaceBoundaryTypes->isDirichlet(Indices::velocity(directionIndex)))
                return problem.dirichlet(element, staggeredLateralScvf)[Indices::velocity(directionIndex)];
            else if (lateralFaceBoundaryTypes->isBeaversJoseph(Indices::velocity(directionIndex)))
            {
                return beaversJosephVelocityAtLateralScvf(problem, element, staggeredFVGeometry, staggeredLateralScvf, faceVars,
                                                          currentScvfBoundaryTypes, lateralFaceBoundaryTypes); // TODO revise BJS signature
            }
            else
                return innerParallelVelocity; // assume a gradient of zero if the boundary conditions don't allow to calculate a gradient
        }();

        // The velocity gradient already accounts for the orientation
        // of the staggered face's outer normal vector. This also correctly accounts for the reduced
        // distance used in the gradient if the lateral scvf lies on a boundary.
        return (outerParallelVelocity - innerParallelVelocity)
               / staggeredFVGeometry.normalDistanceForGradient(staggeredLateralScvf) * staggeredLateralScvf.directionSign();
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
     *                      |  |  |  |  |  ^
     *                      |  |  |  |  |  |
     *
     *              ----------------
     *              |              |
     *              |    in.norm.  |
     *              |       vel.   |
     *              |       ^      |        ^ out.norm.vel.
     *              |       |      |        |
     *       scvf   ---------######O:::::::::       || and # staggered half-control-volume (own element)
     *              |      ||      | curr. ::
     *              |      ||      | scvf  ::       :: staggered half-control-volume (neighbor element)
     *              |      ||      x~~~~>  ::
     *              |      ||      | vel.  ::       # lateral staggered faces over which fluxes are calculated
     *        scvf  |      ||      | Self  ::
     *              ---------#######:::::::::       x dof position
     *                 scvf
     *                                              -- elements
     *
     *                                              O position at which gradient is evaluated
     * \endverbatim
     */
    template<class Problem, class StaggeredFVElementGeometry, class FaceVariables>
    static Scalar velocityGradJI(const Problem& problem,
                                 const Element& element,
                                 const StaggeredFVElementGeometry staggeredFVGeometry,
                                 const typename StaggeredFVElementGeometry::StaggeredSubControlVolumeFace& staggeredLateralScvf,
                                 const FaceVariables& faceVars,
                                 const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
                                 const std::optional<BoundaryTypes>& lateralFaceBoundaryTypes)
    {
        const auto& staggeredScv = staggeredFVGeometry.scv(staggeredLateralScvf.insideScvIdx());
        const auto localLateralFaceIdx = staggeredFVGeometry.localLateralFaceIndex(staggeredLateralScvf); // goes from 0 to dim-1

        // For the velocityGrad_ji gradient, get the velocities perpendicular to the velocity at the current scvf.
        // The inner one is located at staggered face within the own element,
        // the outer one at the respective staggered face of the element on the other side of the
        // current scvf.
        const Scalar innerLateralVelocity = faceVars.velocityLateralInside(localLateralFaceIdx);
        const Scalar outerLateralVelocity = [&]()
        {
            if (!staggeredScv.boundary())
                return faceVars.velocityLateralOutside(localLateralFaceIdx);
            else if (currentScvfBoundaryTypes->isDirichlet(Indices::velocity(staggeredLateralScvf.directionIndex())))
                return problem.dirichlet(element, staggeredLateralScvf)[Indices::velocity(staggeredLateralScvf.directionIndex())];
            else if (currentScvfBoundaryTypes->isBeaversJoseph(Indices::velocity(staggeredLateralScvf.directionIndex())))
            {
                return beaversJosephVelocityAtCurrentScvf(problem, element, staggeredFVGeometry, staggeredLateralScvf, faceVars,
                                                          currentScvfBoundaryTypes, lateralFaceBoundaryTypes);
            }
            else
                return innerLateralVelocity; // assume a gradient of zero if the boundary conditions don't allow to calculate a gradient
        }();

        // Calculate the velocity gradient in positive coordinate direction.
        const Scalar lateralDeltaV = !std::signbit(staggeredScv.directionSign())
                                    ? (outerLateralVelocity - innerLateralVelocity)
                                    : (innerLateralVelocity - outerLateralVelocity);

        return lateralDeltaV / staggeredFVGeometry.tangentialDistanceForGradient(staggeredLateralScvf); // TODO put sign in distance
    }

    /*!
     * \brief Returns the Beavers-Jospeh slip velocity for a scvf which lies on the boundary itself.
     *
     * \verbatim
     *                  in.norm.  B-J slip
     *                     vel.   vel.
     *                     ^       ^
     *                     |       |
     *       scvf   ---------######|*               * boundary
     *              |      ||      |* curr.
     *              |      ||      |* scvf          || and # staggered half-control-volume (own element)
     *              |      ||      x~~~~>
     *              |      ||      |* vel.          # lateral staggered faces
     *        scvf  |      ||      |* Self
     *              ---------#######*                x dof position
     *                 scvf
     *                                              -- element
     * \endverbatim
     *
     */
    template<class Problem, class StaggeredFVElementGeometry, class FaceVariables>
    static Scalar beaversJosephVelocityAtCurrentScvf(const Problem& problem,
                                                     const Element& element,
                                                     const StaggeredFVElementGeometry staggeredFVGeometry,
                                                     const typename StaggeredFVElementGeometry::StaggeredSubControlVolumeFace& staggeredScvf,
                                                     const FaceVariables& faceVars,
                                                     const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
                                                     const std::optional<BoundaryTypes>& lateralFaceBoundaryTypes)
    {
        assert(staggeredScvf.isLateral());
        const auto& staggeredScv = staggeredFVGeometry.scv(staggeredScvf.insideScvIdx());

        const auto tangentialVelocityGradient = [&]()
        {
            // If the current scvf is on a boundary and if a Dirichlet BC for the pressure or a BJ condition for
            // the slip velocity is set there, assume a tangential velocity gradient of zero along the lateral face
            // (towards the current scvf).
            static const bool unsymmetrizedGradientForBJ = getParamFromGroup<bool>(problem.paramGroup(),
                                                           "FreeFlow.EnableUnsymmetrizedVelocityGradientForBeaversJoseph", false);

            if (unsymmetrizedGradientForBJ)
                return 0.0;
            else if (staggeredScvf.boundary() && lateralFaceBoundaryTypes->isBeaversJoseph(Indices::velocity(staggeredScv.directionIndex())))
                return 0.0;
            else
                return velocityGradIJ(problem, element, staggeredFVGeometry, staggeredScvf, faceVars, currentScvfBoundaryTypes, lateralFaceBoundaryTypes);
        }();

        for (const auto& scvfOnBoundary : scvfs(staggeredFVGeometry))
        {
            if (scvfOnBoundary.isFrontal() && scvfOnBoundary.boundary())
            {
                GlobalPosition orientation(0.0);
                orientation[staggeredScvf.directionIndex()] = 1.0;
                const Scalar innerLateralVelocity = faceVars.velocityLateralInside(staggeredFVGeometry.localLateralFaceIndex(staggeredScvf));
                const Scalar distanceNormalToBoundary = staggeredFVGeometry.tangentialDistanceForGradient(staggeredScvf);

                return problem.beaversJosephVelocity(element,
                                                     scvfOnBoundary,
                                                     orientation,
                                                     innerLateralVelocity,
                                                     distanceNormalToBoundary,
                                                     tangentialVelocityGradient);
            }
        }

        DUNE_THROW(Dune::InvalidStateException, "No boundary scvf found"); // TODO make convenience function in fvGeometry

    }

    /*!
     * \brief Returns the Beavers-Jospeh slip velocity for a lateral scvf which lies on the boundary.
     *
     * \verbatim
     *                             B-J slip                  * boundary
     *              ************** vel. *****
     *       scvf   ---------##### ~~~~> ::::                || and # staggered half-control-volume (own element)
     *              |      ||      | curr. ::
     *              |      ||      | scvf  ::                :: staggered half-control-volume (neighbor element)
     *              |      ||      x~~~~>  ::
     *              |      ||      | vel.  ::                # lateral staggered faces
     *        scvf  |      ||      | Self  ::
     *              ---------#######:::::::::                x dof position
     *                 scvf
     *                                                       -- elements
     * \endverbatim
     */
    template<class Problem, class StaggeredFVElementGeometry, class FaceVariables>
    static Scalar beaversJosephVelocityAtLateralScvf(const Problem& problem,
                                                     const Element& element,
                                                     const StaggeredFVElementGeometry staggeredFVGeometry,
                                                     const typename StaggeredFVElementGeometry::StaggeredSubControlVolumeFace& staggeredScvf,
                                                     const FaceVariables& faceVars,
                                                     const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
                                                     const std::optional<BoundaryTypes>& lateralFaceBoundaryTypes)
    {
        assert(staggeredScvf.isLateral());
        const auto& staggeredScv = staggeredFVGeometry.scv(staggeredScvf.insideScvIdx());

        const auto tangentialVelocityGradient = [&]()
        {
            // If the current scvf is on a boundary and if a Dirichlet BC for the pressure or a BJ condition for
            // the slip velocity is set there, assume a tangential velocity gradient of zero along the lateral face
            // (towards the current scvf).
            static const bool unsymmetrizedGradientForBJ = getParamFromGroup<bool>(problem.paramGroup(),
                                                           "FreeFlow.EnableUnsymmetrizedVelocityGradientForBeaversJoseph", false);

            if (unsymmetrizedGradientForBJ)
                return 0.0;
            else if (staggeredScv.boundary() && currentScvfBoundaryTypes->isBeaversJoseph(Indices::velocity(staggeredScvf.directionIndex())))
                return 0.0;
            else
                return velocityGradJI(problem, element, staggeredFVGeometry, staggeredScvf, faceVars, currentScvfBoundaryTypes, lateralFaceBoundaryTypes);
        }();

        GlobalPosition orientation(0.0);
        orientation[staggeredScv.directionIndex()] = 1.0;
        const Scalar innerParallelVelocity = faceVars.velocitySelf();
        const Scalar distanceNormalToBoundary = staggeredFVGeometry.normalDistanceForGradient(staggeredScvf);

        return problem.beaversJosephVelocity(element,
                                             staggeredScvf,
                                             orientation,
                                             innerParallelVelocity,
                                             distanceNormalToBoundary,
                                             tangentialVelocityGradient);
    }

};

} // end namespace Dumux

#endif // DUMUX_NAVIERSTOKES_STAGGERED_VELOCITYGRADIENTS_HH
